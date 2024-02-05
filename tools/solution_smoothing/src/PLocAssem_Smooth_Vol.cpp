#include "PLocAssem_Smooth_Vol.hpp"

PLocAssem_Smooth_Vol::PLocAssem_Smooth_Vol(
    IMaterialModel * const &in_matmodel,
    const int &in_isol_dof, const int &in_osol_dof,
    const int &in_nlocbas )
: isol_dof( in_isol_dof ), osol_dof( in_osol_dof ),
  nLocBas( in_nlocbas ),
  vec_size( osol_dof * in_nlocbas )
{
  matmodel = in_matmodel;

  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  Zero_Tangent_Residual();

  print_info();
}

PLocAssem_Smooth_Vol::~PLocAssem_Smooth_Vol()
{
  delete [] Tangent; Tangent = nullptr;
  delete [] Residual; Residual = nullptr;
}

void PLocAssem_Smooth_Vol::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("  Three-dimensional elastodynamics equation: \n");
  if(nLocBas == 4)
    SYS_T::commPrint("  FEM: 4-node Tetrahedral element \n");
  else if(nLocBas == 10)
    SYS_T::commPrint("  FEM: 10-node Tetrahedral element \n");
  else if(nLocBas == 8)
    SYS_T::commPrint("  FEM: 8-node Hexahedral element \n");
  else if(nLocBas == 27)
    SYS_T::commPrint("  FEM: 27-node Hexahedral element \n");
  else SYS_T::print_fatal("Error: unknown elem type.\n");
  SYS_T::print_sep_line();
}

void PLocAssem_Smooth_Vol::Assem_Residual(
    const double * const &isol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad)
{
  const int nqp = quad -> get_num_quadPts();

  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  Zero_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    double ux_x = 0.0, ux_y = 0.0, ux_z = 0.0;
    double uy_x = 0.0, uy_y = 0.0, uy_z = 0.0;
    double uz_x = 0.0, uz_y = 0.0, uz_z = 0.0;

    element->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int iiisol = isol_dof * ii;

      ux_x += isol[iiisol  ] * dR_dx[ii];
      uy_x += isol[iiisol+1] * dR_dx[ii];
      uz_x += isol[iiisol+2] * dR_dx[ii];

      ux_y += isol[iiisol  ] * dR_dy[ii];
      uy_y += isol[iiisol+1] * dR_dy[ii];
      uz_y += isol[iiisol+2] * dR_dy[ii];

      ux_z += isol[iiisol  ] * dR_dz[ii];
      uy_z += isol[iiisol+1] * dR_dz[ii];
      uz_z += isol[iiisol+2] * dR_dz[ii];
    }

    const Tensor2_3D F(ux_x, ux_y, ux_z,
                       uy_x, uy_y, uy_z,
                       uz_x, uz_y, uz_z);
    
    const Tensor2_3D stress = matmodel->get_Cauchy_stress( F );

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    
    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A];

      Residual[osol_dof*A  ] += gwts * NA * stress(0);
      Residual[osol_dof*A+1] += gwts * NA * stress(4);
      Residual[osol_dof*A+2] += gwts * NA * stress(8);
      Residual[osol_dof*A+3] += gwts * NA * stress(1);
      Residual[osol_dof*A+4] += gwts * NA * stress(5);
      Residual[osol_dof*A+5] += gwts * NA * stress(2);
    }
  }
}

void PLocAssem_Smooth_Vol::Assem_Mass_Residual(
    const double * const &isol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  const int nqp = quad -> get_num_quadPts();

  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  Zero_Tangent_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    double ux_x = 0.0, ux_y = 0.0, ux_z = 0.0;
    double uy_x = 0.0, uy_y = 0.0, uy_z = 0.0;
    double uz_x = 0.0, uz_y = 0.0, uz_z = 0.0;

    element->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int iiisol = isol_dof * ii;

      ux_x += isol[iiisol  ] * dR_dx[ii];
      uy_x += isol[iiisol+1] * dR_dx[ii];
      uz_x += isol[iiisol+2] * dR_dx[ii];

      ux_y += isol[iiisol  ] * dR_dy[ii];
      uy_y += isol[iiisol+1] * dR_dy[ii];
      uz_y += isol[iiisol+2] * dR_dy[ii];

      ux_z += isol[iiisol  ] * dR_dz[ii];
      uy_z += isol[iiisol+1] * dR_dz[ii];
      uz_z += isol[iiisol+2] * dR_dz[ii];
    }

    const Tensor2_3D F(ux_x, ux_y, ux_z,
                       uy_x, uy_y, uy_z,
                       uz_x, uz_y, uz_z);
    
    const Tensor2_3D stress = matmodel->get_Cauchy_stress( F );

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    
    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A];

      Residual[osol_dof*A  ] += gwts * NA * stress(0);
      Residual[osol_dof*A+1] += gwts * NA * stress(4);
      Residual[osol_dof*A+2] += gwts * NA * stress(8);
      Residual[osol_dof*A+3] += gwts * NA * stress(1);
      Residual[osol_dof*A+4] += gwts * NA * stress(5);
      Residual[osol_dof*A+5] += gwts * NA * stress(2);

      for(int B=0; B<nLocBas; ++B)
      {
        // Non-zero diagonal of each A-B block
        for(int C=0; C<osol_dof; ++C)
        {
          Tangent[osol_dof*nLocBas*(osol_dof*A + C) + osol_dof*B + C] += gwts * NA * R[B];
        }
      }
    }
  }
}

//EOF