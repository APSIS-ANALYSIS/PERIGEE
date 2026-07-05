#include "PGAssem_Transport_NS_FEM.hpp"

PGAssem_Transport_NS_FEM::PGAssem_Transport_NS_FEM(
    const ALocal_InflowBC * const in_inflow,
    std::unique_ptr<ALocal_EBC> in_ebc,
    std::unique_ptr<ALocal_IEN> in_locien,
    std::unique_ptr<ALocal_Elem> in_locelem,
    std::unique_ptr<FEANode> in_fnode,
    std::unique_ptr<APart_Node> in_pnode,
    const FEType &in_elemType,
    const int &in_nqp_s,
    std::unique_ptr<PLocAssem_Transport_VMS_NS_GenAlpha> in_locassem,
    const int &in_nz_estimate )
: inflow(in_inflow),
  ebc(std::move(in_ebc)),
  locien(std::move(in_locien)),
  locelem(std::move(in_locelem)),
  fnode(std::move(in_fnode)),
  pnode(std::move(in_pnode)),
  elemType(in_elemType),
  nqps(in_nqp_s),
  elements(ElementFactory::createSurElement(elemType, nqps)),
  quads(QuadPtsFactory::createSurQuadrature(elemType, nqps)),
  locassem(std::move(in_locassem)),
  nLocBas(this->locassem->get_nLocBas()),
  snLocBas(elements->get_nLocBas()),
  dof_mat(this->locassem->get_dof()),
  nlgn(this->pnode->get_nlocghonode())
{
  const int nlocrow = pnode->get_nlocalnode() * dof_mat;

  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat * in_nz_estimate, NULL,
      dof_mat * in_nz_estimate, NULL, &K);

  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, nlocrow, PETSC_DECIDE);
  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  Release_nonzero_err_str();

  Assem_nonzero_estimate();

  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);
  MatDestroy(&K);

  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}

PGAssem_Transport_NS_FEM::~PGAssem_Transport_NS_FEM()
{
  VecDestroy(&G);
  MatDestroy(&K);
}

void PGAssem_Transport_NS_FEM::Apply_inflow_BC_G()
{
  for(int face=0; face<inflow->get_num_nbc(); ++face)
  {
    const int nnode = inflow->get_Num_LD(face);
    for(int ii=0; ii<nnode; ++ii)
      VecSetValue(G, inflow->get_LDN(face, ii), 0.0, INSERT_VALUES);
  }
}

void PGAssem_Transport_NS_FEM::Apply_inflow_BC_KG()
{
  for(int face=0; face<inflow->get_num_nbc(); ++face)
  {
    const int nnode = inflow->get_Num_LD(face);
    for(int ii=0; ii<nnode; ++ii)
    {
      const int row = inflow->get_LDN(face, ii);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, row, 1.0, ADD_VALUES);
    }
  }
}

void PGAssem_Transport_NS_FEM::Assem_nonzero_estimate()
{
  const int nElem = locelem->get_nlocalele();
  locassem->Assem_Estimate();

  PetscInt * row_index = new PetscInt[nLocBas];
  for(int ee=0; ee<nElem; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii)
      row_index[ii] = pnode->get_local_to_global( locien->get_LIEN(ee, ii) );

    MatSetValues(K, nLocBas, row_index, nLocBas, row_index, locassem->Tangent, ADD_VALUES);
  }

  delete [] row_index; row_index = nullptr;

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
  Apply_inflow_BC_KG();
  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

double PGAssem_Transport_NS_FEM::Assem_surface_HI_flowrate(
    const PDNSolution * const &transport_sol,
    const PDNSolution * const &flow_sol,
    const int &ebc_id ) const
{
  SYS_T::print_fatal_if(!has_hi_model(),
      "Error: Assem_surface_HI_flowrate requires a hemolysis model.\n");

  double * transport_array = new double[nlgn];
  double * flow_array = new double[nlgn * 4];
  double * local_transport = new double[snLocBas];
  double * local_flow = new double[snLocBas * 4];
  int * LSIEN = new int[snLocBas];
  double * sctrl_x = new double[snLocBas];
  double * sctrl_y = new double[snLocBas];
  double * sctrl_z = new double[snLocBas];

  transport_sol->GetLocalArray(transport_array);
  flow_sol->GetLocalArray(flow_array);

  const int num_sele = ebc->get_num_local_cell(ebc_id);
  const double beta = get_hi_beta();
  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc->get_SIEN(ebc_id, ee, LSIEN);
    ebc->get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      local_transport[ii] = transport_array[LSIEN[ii]];

      for(int jj=0; jj<4; ++jj)
        local_flow[4*ii + jj] = flow_array[4*LSIEN[ii] + jj];
    }

    elements->buildBasis(quads.get(), sctrl_x, sctrl_y, sctrl_z);

    for(int qua=0; qua<nqps; ++qua)
    {
      const std::vector<double> R = elements->get_R(qua);
      double surface_area;
      const Vector_3 n_out = elements->get_2d_normal_out(qua, surface_area);

      double D = 0.0;
      Vector_3 velo(0.0, 0.0, 0.0);
      for(int ii=0; ii<snLocBas; ++ii)
      {
        D += local_transport[ii] * R[ii];
        velo.x() += local_flow[4*ii+1] * R[ii];
        velo.y() += local_flow[4*ii+2] * R[ii];
        velo.z() += local_flow[4*ii+3] * R[ii];
      }

      const double HI = std::pow(std::max(0.0, D), beta);
      esum += surface_area * quads->get_qw(qua) * HI * Vec3::dot_product(velo, n_out);
    }
  }

  delete [] transport_array;
  delete [] flow_array;
  delete [] local_transport;
  delete [] local_flow;
  delete [] LSIEN;
  delete [] sctrl_x;
  delete [] sctrl_y;
  delete [] sctrl_z;

  double sum = 0.0;
  MPI_Allreduce(&esum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  return sum;
}

double PGAssem_Transport_NS_FEM::Assem_surface_HI_integral(
    const PDNSolution * const &transport_sol,
    const int &ebc_id ) const
{
  SYS_T::print_fatal_if(!has_hi_model(),
      "Error: Assem_surface_HI_integral requires a hemolysis model.\n");

  double * transport_array = new double[nlgn];
  double * local_transport = new double[snLocBas];
  int * LSIEN = new int[snLocBas];
  double * sctrl_x = new double[snLocBas];
  double * sctrl_y = new double[snLocBas];
  double * sctrl_z = new double[snLocBas];

  transport_sol->GetLocalArray(transport_array);

  const int num_sele = ebc->get_num_local_cell(ebc_id);
  const double beta = get_hi_beta();
  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc->get_SIEN(ebc_id, ee, LSIEN);
    ebc->get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    for(int ii=0; ii<snLocBas; ++ii) local_transport[ii] = transport_array[LSIEN[ii]];

    elements->buildBasis(quads.get(), sctrl_x, sctrl_y, sctrl_z);

    for(int qua=0; qua<nqps; ++qua)
    {
      const std::vector<double> R = elements->get_R(qua);
      double surface_area;
      elements->get_2d_normal_out(qua, surface_area);

      double D = 0.0;
      for(int ii=0; ii<snLocBas; ++ii) D += local_transport[ii] * R[ii];

      const double HI = std::pow(std::max(0.0, D), beta);
      esum += surface_area * quads->get_qw(qua) * HI;
    }
  }

  delete [] transport_array;
  delete [] local_transport;
  delete [] LSIEN;
  delete [] sctrl_x;
  delete [] sctrl_y;
  delete [] sctrl_z;

  double sum = 0.0;
  MPI_Allreduce(&esum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  return sum;
}

void PGAssem_Transport_NS_FEM::Assem_mass_residual(
    const PDNSolution * const &phi,
    const PDNSolution * const &flow_sol )
{
  const int nElem = locelem->get_nlocalele();

  double * phi_array = new double[nlgn];
  double * flow_array = new double[nlgn * 4];
  double * local_phi = new double[nLocBas];
  double * local_flow = new double[nLocBas * 4];
  int * IEN_e = new int[nLocBas];
  double * ectrl_x = new double[nLocBas];
  double * ectrl_y = new double[nLocBas];
  double * ectrl_z = new double[nLocBas];
  PetscInt * row_index = new PetscInt[nLocBas];

  phi->GetLocalArray(phi_array);
  flow_sol->GetLocalArray(flow_array);

  for(int ee=0; ee<nElem; ++ee)
  {
    locien->get_LIEN(ee, IEN_e);
    GetLocal(phi_array, IEN_e, 1, local_phi);
    GetLocal(flow_array, IEN_e, 4, local_flow);
    fnode->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    locassem->Assem_Mass_Residual(local_phi, local_flow, ectrl_x, ectrl_y, ectrl_z);

    for(int ii=0; ii<nLocBas; ++ii)
      row_index[ii] = pnode->get_local_to_global(IEN_e[ii]);

    MatSetValues(K, nLocBas, row_index, nLocBas, row_index, locassem->Tangent, ADD_VALUES);
    VecSetValues(G, nLocBas, row_index, locassem->Residual, ADD_VALUES);
  }

  delete [] phi_array;
  delete [] flow_array;
  delete [] local_phi;
  delete [] local_flow;
  delete [] IEN_e;
  delete [] ectrl_x;
  delete [] ectrl_y;
  delete [] ectrl_z;
  delete [] row_index;

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
  Apply_inflow_BC_KG();
  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_Transport_NS_FEM::Assem_residual(
    const PDNSolution * const &dot_phi,
    const PDNSolution * const &phi,
    const PDNSolution * const &flow_sol,
    const double &curr_time,
    const double &dt )
{
  const int nElem = locelem->get_nlocalele();

  double * dot_array = new double[nlgn];
  double * phi_array = new double[nlgn];
  double * flow_array = new double[nlgn * 4];
  double * local_dot = new double[nLocBas];
  double * local_phi = new double[nLocBas];
  double * local_flow = new double[nLocBas * 4];
  int * IEN_e = new int[nLocBas];
  double * ectrl_x = new double[nLocBas];
  double * ectrl_y = new double[nLocBas];
  double * ectrl_z = new double[nLocBas];
  PetscInt * row_index = new PetscInt[nLocBas];

  dot_phi->GetLocalArray(dot_array);
  phi->GetLocalArray(phi_array);
  flow_sol->GetLocalArray(flow_array);

  for(int ee=0; ee<nElem; ++ee)
  {
    locien->get_LIEN(ee, IEN_e);
    GetLocal(dot_array, IEN_e, 1, local_dot);
    GetLocal(phi_array, IEN_e, 1, local_phi);
    GetLocal(flow_array, IEN_e, 4, local_flow);
    fnode->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    locassem->Assem_Residual(curr_time, dt, local_dot, local_phi, local_flow,
        ectrl_x, ectrl_y, ectrl_z);

    for(int ii=0; ii<nLocBas; ++ii)
      row_index[ii] = pnode->get_local_to_global(IEN_e[ii]);

    VecSetValues(G, nLocBas, row_index, locassem->Residual, ADD_VALUES);
  }

  delete [] dot_array;
  delete [] phi_array;
  delete [] flow_array;
  delete [] local_dot;
  delete [] local_phi;
  delete [] local_flow;
  delete [] IEN_e;
  delete [] ectrl_x;
  delete [] ectrl_y;
  delete [] ectrl_z;
  delete [] row_index;

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
  Apply_inflow_BC_G();
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_Transport_NS_FEM::Assem_tangent_residual(
    const PDNSolution * const &dot_phi,
    const PDNSolution * const &phi,
    const PDNSolution * const &flow_sol,
    const double &curr_time,
    const double &dt )
{
  const int nElem = locelem->get_nlocalele();

  double * dot_array = new double[nlgn];
  double * phi_array = new double[nlgn];
  double * flow_array = new double[nlgn * 4];
  double * local_dot = new double[nLocBas];
  double * local_phi = new double[nLocBas];
  double * local_flow = new double[nLocBas * 4];
  int * IEN_e = new int[nLocBas];
  double * ectrl_x = new double[nLocBas];
  double * ectrl_y = new double[nLocBas];
  double * ectrl_z = new double[nLocBas];
  PetscInt * row_index = new PetscInt[nLocBas];

  dot_phi->GetLocalArray(dot_array);
  phi->GetLocalArray(phi_array);
  flow_sol->GetLocalArray(flow_array);

  for(int ee=0; ee<nElem; ++ee)
  {
    locien->get_LIEN(ee, IEN_e);
    GetLocal(dot_array, IEN_e, 1, local_dot);
    GetLocal(phi_array, IEN_e, 1, local_phi);
    GetLocal(flow_array, IEN_e, 4, local_flow);
    fnode->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    locassem->Assem_Tangent_Residual(curr_time, dt, local_dot, local_phi,
        local_flow, ectrl_x, ectrl_y, ectrl_z);

    for(int ii=0; ii<nLocBas; ++ii)
      row_index[ii] = pnode->get_local_to_global(IEN_e[ii]);

    MatSetValues(K, nLocBas, row_index, nLocBas, row_index, locassem->Tangent, ADD_VALUES);
    VecSetValues(G, nLocBas, row_index, locassem->Residual, ADD_VALUES);
  }

  delete [] dot_array;
  delete [] phi_array;
  delete [] flow_array;
  delete [] local_dot;
  delete [] local_phi;
  delete [] local_flow;
  delete [] IEN_e;
  delete [] ectrl_x;
  delete [] ectrl_y;
  delete [] ectrl_z;
  delete [] row_index;

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
  Apply_inflow_BC_KG();
  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}
