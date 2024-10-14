#include "Sliding_Interface_Tools.hpp"

namespace SI_T
{
  SI_solution::SI_solution(const std::string &fileBaseName, const int &in_cpu_rank)
    : cpu_rank( in_cpu_rank )
  {
    const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

    hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

    HDF5_Reader * h5r = new HDF5_Reader( file_id );

    const std::string gname("/sliding");

    const int num_itf = h5r -> read_intScalar( gname.c_str(), "num_interface" );

    num_fixed_node.assign(num_itf, 0);
    fixed_node_sol.resize(num_itf);

    f_r_node_sol.resize(num_itf);
    f_r_node_sol_x.resize(num_itf);
    f_r_node_sol_y.resize(num_itf);
    f_r_node_sol_z.resize(num_itf);

    f_r_node_mvelo.resize(num_itf);

    L2_proj_fixed_node_pos.resize(num_itf);
    num_all_fixed_node = h5r -> read_intScalar( gname.c_str(), "num_all_fixed_node" );
    all_fixed_inner_node = h5r -> read_intVector( gname.c_str(), "all_fixed_inner_node" );
    num_all_fixed_inner_node = VEC_T::get_size(all_fixed_inner_node);

    fixed_node_part_tag.resize(num_itf);
    fixed_node_loc_pos.resize(num_itf);

    num_rotated_node.assign(num_itf, 0);
    rotated_node_sol.resize(num_itf);

    r_f_node_sol.resize(num_itf);
    r_f_node_sol_x.resize(num_itf);
    r_f_node_sol_y.resize(num_itf);
    r_f_node_sol_z.resize(num_itf);

    L2_proj_rotated_node_pos.resize(num_itf);
    num_all_rotated_node = h5r -> read_intScalar( gname.c_str(), "num_all_rotated_node" );
    all_rotated_inner_node = h5r -> read_intVector( gname.c_str(), "all_rotated_inner_node" );
    num_all_rotated_inner_node = VEC_T::get_size(all_rotated_inner_node);

    rotated_node_mvelo.resize(num_itf);
    rotated_node_mdisp.resize(num_itf);
    rotated_node_part_tag.resize(num_itf);
    rotated_node_loc_pos.resize(num_itf);

    const std::string mesh_info("/Global_Mesh_Info");
    nLocBas = h5r -> read_intScalar(mesh_info.c_str(), "nLocBas");
    dof_sol = h5r -> read_intScalar(mesh_info.c_str(), "dofNum");

    std::string groupbase(gname);
    groupbase.append("/interfaceid_");

    for(int ii=0; ii<num_itf; ++ii)
    {
      std::string subgroup_name(groupbase);
      subgroup_name.append( std::to_string(ii) );
      
      num_fixed_node[ii] = VEC_T::get_size(h5r -> read_intVector( subgroup_name.c_str(), "fixed_node_map" ));

      fixed_node_sol[ii] = std::vector<double> (dof_sol * num_fixed_node[ii], 0.0);

      f_r_node_sol[ii] = std::vector<double> (dof_sol * num_fixed_node[ii], 0.0);

      f_r_node_sol_x[ii] = std::vector<double> (3 * num_fixed_node[ii], 0.0);
      f_r_node_sol_y[ii] = std::vector<double> (3 * num_fixed_node[ii], 0.0);
      f_r_node_sol_z[ii] = std::vector<double> (3 * num_fixed_node[ii], 0.0);

      f_r_node_mvelo[ii] = std::vector<double> (3 * num_fixed_node[ii], 0.0);

      L2_proj_fixed_node_pos[ii] = h5r -> read_intVector( subgroup_name.c_str(), "L2_proj_fixed_node_pos" );

      fixed_node_part_tag[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_node_part_tag" );

      fixed_node_loc_pos[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_node_loc_pos" );

      num_rotated_node[ii] = VEC_T::get_size(h5r -> read_intVector( subgroup_name.c_str(), "rotated_node_map" ));

      rotated_node_sol[ii] = std::vector<double> (dof_sol * num_rotated_node[ii], 0.0);

      r_f_node_sol[ii] = std::vector<double> (dof_sol * num_rotated_node[ii], 0.0);

      r_f_node_sol_x[ii] = std::vector<double> (3 * num_rotated_node[ii], 0.0);
      r_f_node_sol_y[ii] = std::vector<double> (3 * num_rotated_node[ii], 0.0);
      r_f_node_sol_z[ii] = std::vector<double> (3 * num_rotated_node[ii], 0.0);

      L2_proj_rotated_node_pos[ii] = h5r -> read_intVector( subgroup_name.c_str(), "L2_proj_rotated_node_pos" );

      rotated_node_mvelo[ii] = std::vector<double> (dof_sol * num_rotated_node[ii], 0.0);

      rotated_node_mdisp[ii] = std::vector<double> (dof_sol * num_rotated_node[ii], 0.0); 

      rotated_node_part_tag[ii] = h5r -> read_intVector( subgroup_name.c_str(), "rotated_node_part_tag" );

      rotated_node_loc_pos[ii] = h5r -> read_intVector( subgroup_name.c_str(), "rotated_node_loc_pos" );
    }

    delete h5r; H5Fclose( file_id );
  }

  void SI_solution::update_node_sol(const PDNSolution * const &sol)
  {
    const int nlgn = sol->get_nlgn();
    double * array = new double [nlgn];

    sol->GetLocalArray( array );

    zero_node_sol();

    for(int ii = 0; ii < VEC_T::get_size(num_fixed_node); ++ii)
    { 
      std::vector<double> temp_fixed_node_sol = fixed_node_sol[ii];
      for(int nn = 0; nn < num_fixed_node[ii]; ++nn)
      {
        // Pick out the solution of nodes in each part
        if(fixed_node_part_tag[ii][nn] == cpu_rank)
        {  
          // Like GetLocal in PGAssem
          const int loc_pos = fixed_node_loc_pos[ii][nn];
          for(int dd = 0; dd < dof_sol; ++dd)
            temp_fixed_node_sol[dof_sol * nn + dd] = array[dof_sol * loc_pos + dd];
        }
      }

      std::vector<double> temp_rotated_node_sol = rotated_node_sol[ii];
      for(int nn = 0; nn < num_rotated_node[ii]; ++nn)
      {
        // Pick out the solution of nodes in each part
        if(rotated_node_part_tag[ii][nn] == cpu_rank)
        { 
          // Like GetLocal in PGAssem
          const int loc_pos = rotated_node_loc_pos[ii][nn];
          for(int dd = 0; dd < dof_sol; ++dd)
            temp_rotated_node_sol[dof_sol * nn + dd] = array[dof_sol * loc_pos + dd];
        }
      }

      // Summation from each part
      MPI_Allreduce(&temp_fixed_node_sol[0], &fixed_node_sol[ii][0], dof_sol * num_fixed_node[ii], MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
      MPI_Allreduce(&temp_rotated_node_sol[0], &rotated_node_sol[ii][0], dof_sol * num_rotated_node[ii], MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    }

    delete [] array; array = nullptr;
  }

  void SI_solution::update_node_mvelo(const PDNSolution * const &mvelo)
  {
    const int nlgn = mvelo->get_nlgn();    
    double * array = new double [nlgn];

    mvelo->GetLocalArray( array );

    Zero_node_mvelo();

    for(int ii = 0; ii < VEC_T::get_size(num_fixed_node); ++ii)
    { 
      std::vector<double> temp_rotated_node_mvelo = rotated_node_mvelo[ii];
      for(int nn = 0; nn < num_rotated_node[ii]; ++nn)
      {
        // Pick out the solution of nodes in each part
        if(rotated_node_part_tag[ii][nn] == cpu_rank)
        { 
          // Like GetLocal in PGAssem
          const int loc_pos = rotated_node_loc_pos[ii][nn];
          for(int dd = 0; dd < 3; ++dd)
            temp_rotated_node_mvelo[3 * nn + dd] = array[3 * loc_pos + dd];
        }
      }

      // Summation from each part
      MPI_Allreduce(&temp_rotated_node_mvelo[0], &rotated_node_mvelo[ii][0], 3 * num_rotated_node[ii], MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    }

    delete [] array; array = nullptr;
  }

  void SI_solution::update_node_mdisp(const PDNSolution * const &mdisp)
  {
    const int nlgn = mdisp->get_nlgn();    
    double * array = new double [nlgn];

    mdisp->GetLocalArray( array );

    Zero_node_disp();

    for(int ii = 0; ii < VEC_T::get_size(num_fixed_node); ++ii)
    { 
      std::vector<double> temp_rotated_node_mdisp = rotated_node_mdisp[ii];
      for(int nn = 0; nn < num_rotated_node[ii]; ++nn)
      {
        // Pick out the solution of nodes in each part
        if(rotated_node_part_tag[ii][nn] == cpu_rank)
        { 
          // Like GetLocal in PGAssem
          const int loc_pos = rotated_node_loc_pos[ii][nn];
          for(int dd = 0; dd < 3; ++dd)
            temp_rotated_node_mdisp[3 * nn + dd] = array[3 * loc_pos + dd];
        }
      }

      // Summation from each part
      MPI_Allreduce(&temp_rotated_node_mdisp[0], &rotated_node_mdisp[ii][0], 3 * num_rotated_node[ii], MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    }

    delete [] array; array = nullptr; 
  }

  void SI_solution::update_L2_proj_result(const Vec &L2_proj_lhs, const int &which)
  {
    for(int ii = 0; ii < VEC_T::get_size(num_fixed_node); ++ii)
    {
      PetscInt * idx;
      if(which == 0)
      {
        idx = new PetscInt[4 * num_fixed_node[ii]];
        for(int jj = 0; jj < num_fixed_node[ii]; ++jj)
        {
          idx[4 * jj] = 4 * L2_proj_fixed_node_pos[ii][jj];
          idx[4 * jj + 1] = 4 * L2_proj_fixed_node_pos[ii][jj] + 1;
          idx[4 * jj + 2] = 4 * L2_proj_fixed_node_pos[ii][jj] + 2;
          idx[4 * jj + 3] = 4 * L2_proj_fixed_node_pos[ii][jj] + 3;
        }
      }
      else if(which > 0 && which < 5)
      {
        idx = new PetscInt[3 * num_fixed_node[ii]];
        for(int jj = 0; jj < num_fixed_node[ii]; ++jj)
        {
          idx[3 * jj] = 4 * L2_proj_fixed_node_pos[ii][jj] + 1;
          idx[3 * jj + 1] = 4 * L2_proj_fixed_node_pos[ii][jj] + 2;
          idx[3 * jj + 2] = 4 * L2_proj_fixed_node_pos[ii][jj] + 3;
        }
      }
      else if(which == 5)
      {
        idx = new PetscInt[4 * num_rotated_node[ii]];
        for(int jj = 0; jj < num_rotated_node[ii]; ++jj)
        {
          idx[4 * jj] = 4 * L2_proj_rotated_node_pos[ii][jj];
          idx[4 * jj + 1] = 4 * L2_proj_rotated_node_pos[ii][jj] + 1;
          idx[4 * jj + 2] = 4 * L2_proj_rotated_node_pos[ii][jj] + 2;
          idx[4 * jj + 3] = 4 * L2_proj_rotated_node_pos[ii][jj] + 3;
        }
      }
      else if(which > 5 && which < 9)
      {
        idx = new PetscInt[3 * num_rotated_node[ii]];
        for(int jj = 0; jj < num_rotated_node[ii]; ++jj)
        {
          idx[3 * jj] = 4 * L2_proj_rotated_node_pos[ii][jj] + 1;
          idx[3 * jj + 1] = 4 * L2_proj_rotated_node_pos[ii][jj] + 2;
          idx[3 * jj + 2] = 4 * L2_proj_rotated_node_pos[ii][jj] + 3;
        }
      }

      switch (which)
      {
        case 0:
          PETSc_T::Scatter(L2_proj_lhs, idx, 4 * num_fixed_node[ii], &f_r_node_sol[ii][0]);
        break;

        case 1:
          PETSc_T::Scatter(L2_proj_lhs, idx, 3 * num_fixed_node[ii], &f_r_node_sol_x[ii][0]);
        break;

        case 2:
          PETSc_T::Scatter(L2_proj_lhs, idx, 3 * num_fixed_node[ii], &f_r_node_sol_y[ii][0]);
        break;

        case 3:
          PETSc_T::Scatter(L2_proj_lhs, idx, 3 * num_fixed_node[ii], &f_r_node_sol_z[ii][0]);
        break;

        case 4:
          PETSc_T::Scatter(L2_proj_lhs, idx, 3 * num_fixed_node[ii], &f_r_node_mvelo[ii][0]);
        break;

        case 5:
          PETSc_T::Scatter(L2_proj_lhs, idx, 4 * num_rotated_node[ii], &r_f_node_sol[ii][0]);
        break;

        case 6:
          PETSc_T::Scatter(L2_proj_lhs, idx, 3 * num_rotated_node[ii], &r_f_node_sol_x[ii][0]);
        break;

        case 7:
          PETSc_T::Scatter(L2_proj_lhs, idx, 3 * num_rotated_node[ii], &r_f_node_sol_y[ii][0]);
        break;

        case 8:
          PETSc_T::Scatter(L2_proj_lhs, idx, 3 * num_rotated_node[ii], &r_f_node_sol_z[ii][0]);
        break;
        
        default:
          SYS_T::print_fatal("Error, SI_solution::update_L2_proj_result, wrong input.\n");
        break;
      }
        
      delete idx; idx = nullptr;
    }
  }

  SI_quad_point::SI_quad_point(const ALocal_Interface * const &itf, const int &nqp_sur_in)
    : nqp_sur( nqp_sur_in )
  {
    // fixed_qp_curr_rotated_tag.resize(itf->get_num_itf());
    fixed_qp_curr_rotated_ee.resize(itf->get_num_itf());
    fixed_qp_curr_rotated_xi.resize(itf->get_num_itf());
    fixed_qp_curr_rotated_eta.resize(itf->get_num_itf());

    // rotated_qp_curr_fixed_tag.resize(itf->get_num_itf());
    rotated_qp_curr_fixed_ee.resize(itf->get_num_itf());
    rotated_qp_curr_fixed_xi.resize(itf->get_num_itf());
    rotated_qp_curr_fixed_eta.resize(itf->get_num_itf());

    for(int ii = 0; ii < itf->get_num_itf(); ++ii)
    {
      // fixed_qp_curr_rotated_tag[ii].assign(itf->get_num_fixed_ele(ii) * nqp_sur, -1);
      fixed_qp_curr_rotated_ee[ii].assign(itf->get_num_fixed_ele(ii) * nqp_sur, -1);
      fixed_qp_curr_rotated_xi[ii].assign(itf->get_num_fixed_ele(ii) * nqp_sur, 0.0);
      fixed_qp_curr_rotated_eta[ii].assign(itf->get_num_fixed_ele(ii) * nqp_sur, 0.0);

      // rotated_qp_curr_fixed_tag[ii].assign(itf->get_num_rotated_ele(ii) * nqp_sur, -1);
      rotated_qp_curr_fixed_ee[ii].assign(itf->get_num_rotated_ele(ii) * nqp_sur, -1);
      rotated_qp_curr_fixed_xi[ii].assign(itf->get_num_rotated_ele(ii) * nqp_sur, 0.0);
      rotated_qp_curr_fixed_eta[ii].assign(itf->get_num_rotated_ele(ii) * nqp_sur, 0.0);
    }
  }

  void SI_quad_point::set_curr_rotated(const double &itf_id, const int &fixed_ee_index, const int &qua,
    const int &rotated_ee, const double &xi, const double &eta)
  {
    // fixed_qp_curr_rotated_tag[itf_id][fixed_ee * nqp_sur + qua]   = ele_tag;
    fixed_qp_curr_rotated_ee[itf_id][fixed_ee_index * nqp_sur + qua]    = rotated_ee;
    fixed_qp_curr_rotated_xi[itf_id][fixed_ee_index * nqp_sur + qua] = xi;
    fixed_qp_curr_rotated_eta[itf_id][fixed_ee_index * nqp_sur + qua] = eta;
  }

  void SI_quad_point::get_curr_rotated(const double &itf_id, const int &fixed_ee_index, const int &qua,
    int &rotated_ee, double &xi, double &eta) const
  {
    // ele_tag = fixed_qp_curr_rotated_tag[itf_id][fixed_ee * nqp_sur + qua];
    rotated_ee = fixed_qp_curr_rotated_ee[itf_id][fixed_ee_index * nqp_sur + qua];
    xi = fixed_qp_curr_rotated_xi[itf_id][fixed_ee_index * nqp_sur + qua];
    eta = fixed_qp_curr_rotated_eta[itf_id][fixed_ee_index * nqp_sur + qua];
  }

  void SI_quad_point::set_curr_fixed(const double &itf_id, const int &rotated_ee_index, const int &qua,
    const int &fixed_ee, const double &xi, const double &eta)
  {
    // rotated_qp_curr_fixed_tag[itf_id][rotated_ee * nqp_sur + qua]   = ele_tag;
    rotated_qp_curr_fixed_ee[itf_id][rotated_ee_index * nqp_sur + qua]    = fixed_ee;
    rotated_qp_curr_fixed_xi[itf_id][rotated_ee_index * nqp_sur + qua] = xi;
    rotated_qp_curr_fixed_eta[itf_id][rotated_ee_index * nqp_sur + qua] = eta;
  }

  void SI_quad_point::get_curr_fixed(const double &itf_id, const int &rotated_ee_index, const int &qua,
    int &fixed_ee, double &xi, double &eta) const
  {
    // ele_tag = rotated_qp_curr_fixed_tag[itf_id][rotated_ee * nqp_sur + qua];
    fixed_ee = rotated_qp_curr_fixed_ee[itf_id][rotated_ee_index * nqp_sur + qua];
    xi = rotated_qp_curr_fixed_xi[itf_id][rotated_ee_index * nqp_sur + qua];
    eta = rotated_qp_curr_fixed_eta[itf_id][rotated_ee_index * nqp_sur + qua];
  }

  void SI_quad_point::search_all_opposite_point(
    FEAElement * const &anchor_elementv,
    FEAElement * const &opposite_elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_s,
    IQuadPts * const &free_quad,
    const ALocal_Interface * const &itf_part,
    const SI_solution * const &SI_sol )
  {
    const int nLocBas = itf_part->get_nLocBas();
    double * ctrl_x = new double [nLocBas];
    double * ctrl_y = new double [nLocBas];
    double * ctrl_z = new double [nLocBas];

    const int num_itf {itf_part->get_num_itf()};

    for(int itf_id{0}; itf_id<num_itf; ++itf_id)
    {
      SYS_T::commPrint("itf_id = %d\n", itf_id);
      const int num_fixed_elem = itf_part->get_num_fixed_ele(itf_id);

      for(int ee_index{0}; ee_index<num_fixed_elem; ++ee_index)
      {
        // SYS_T::commPrint("  fixed_ee = %d\n", ee);
        const int ee = itf_part->get_fixed_ele(itf_id, ee_index);

        itf_part->get_fixed_ele_ctrlPts(itf_id, ee, ctrl_x, ctrl_y, ctrl_z);

        const int fixed_face_id {itf_part->get_fixed_face_id(itf_id, ee)};
        
        anchor_elementv->buildBasis(fixed_face_id, quad_s, ctrl_x, ctrl_y, ctrl_z);

        const int fixed_face_nqp {quad_s->get_num_quadPts()};

        for(int qua{0}; qua<fixed_face_nqp; ++qua)
        {
          auto R = anchor_elementv->get_R(qua);

          // The xyz-coordinates of the quadrature point
          Vector_3 coor(0.0, 0.0, 0.0);
          for(int ii{0}; ii<nLocBas; ++ii)
          {
            coor.x() += ctrl_x[ii] * R[ii];
            coor.y() += ctrl_y[ii] * R[ii];
            coor.z() += ctrl_z[ii] * R[ii];
          }

          // SYS_T::commPrint("    point %d:\n", qua);

          int ele_tag {itf_part->get_fixed_ele_tag(itf_id, ee)};
          int rotated_ee {0};
          search_opposite_rotated_point(coor, itf_part, SI_sol, itf_id, opposite_elementv, elements, ele_tag, rotated_ee, free_quad);

          set_curr_rotated(itf_id, ee_index, qua, rotated_ee, free_quad->get_qp(0, 0), free_quad->get_qp(0, 1));
        }
      }

      const int num_rotated_elem = itf_part->get_num_rotated_ele(itf_id);
      int * rotated_local_ien = new int [nLocBas];
      double * rotated_local_disp = new double [nLocBas * 3];
      double * volctrl_x = new double [nLocBas];
      double * volctrl_y = new double [nLocBas];
      double * volctrl_z = new double [nLocBas];

      for(int ee_index{0}; ee_index<num_rotated_elem; ++ee_index)
      {
        // SYS_T::commPrint("  fixed_ee = %d\n", ee);
        const int ee = itf_part->get_rotated_ele(itf_id, ee_index);

        itf_part->get_rotated_ele_ctrlPts(itf_id, ee, ctrl_x, ctrl_y, ctrl_z);
        SI_sol->get_rotated_mdisp(itf_part, itf_id, ee, rotated_local_ien, rotated_local_disp);
        SI_T::get_currPts(ctrl_x, ctrl_y, ctrl_z, rotated_local_disp, nLocBas, volctrl_x, volctrl_y, volctrl_z);

        const int rotated_face_id {itf_part->get_rotated_face_id(itf_id, ee)};

        anchor_elementv->buildBasis(rotated_face_id, quad_s, volctrl_x, volctrl_y, volctrl_z);

        const int rotated_face_nqp {quad_s->get_num_quadPts()};

        for(int qua{0}; qua<rotated_face_nqp; ++qua)
        {
          auto R = anchor_elementv->get_R(qua);

          // The xyz-coordinates of the quadrature point
          Vector_3 coor(0.0, 0.0, 0.0);
          for(int ii{0}; ii<nLocBas; ++ii)
          {
            coor.x() += volctrl_x[ii] * R[ii];
            coor.y() += volctrl_y[ii] * R[ii];
            coor.z() += volctrl_z[ii] * R[ii];
          }

          // SYS_T::commPrint("    point %d:\n", qua);

          int ele_tag {itf_part->get_rotated_ele_tag(itf_id, ee)};
          int fixed_ee {0};
          search_opposite_fixed_point(coor, itf_part, SI_sol, itf_id, opposite_elementv, elements, ele_tag, fixed_ee, free_quad);

          set_curr_fixed(itf_id, ee_index, qua, fixed_ee, free_quad->get_qp(0, 0), free_quad->get_qp(0, 1));
        }
      }

      delete [] volctrl_x; volctrl_x = nullptr;
      delete [] volctrl_y; volctrl_y = nullptr;
      delete [] volctrl_z; volctrl_z = nullptr;

      delete [] rotated_local_ien; rotated_local_ien = nullptr;
      delete [] rotated_local_disp; rotated_local_disp = nullptr;
    }

    delete [] ctrl_x; ctrl_x = nullptr;
    delete [] ctrl_y; ctrl_y = nullptr;
    delete [] ctrl_z; ctrl_z = nullptr;
  }
  
  void SI_quad_point::search_opposite_rotated_point(
    const Vector_3 &fixed_pt,
    const ALocal_Interface * const &itf_part,
    const SI_solution * const &SI_sol,
    const int &itf_id,
    FEAElement * rotated_elementv,
    FEAElement * elements,
    int &tag,
    int &rotated_ee,
    IQuadPts * const &rotated_xi )
  {
    const int nLocBas = itf_part->get_nLocBas();
    bool is_found = false;

    double * inictrl_x = new double [nLocBas];
    double * inictrl_y = new double [nLocBas];
    double * inictrl_z = new double [nLocBas];

    double * volctrl_x = new double [nLocBas];
    double * volctrl_y = new double [nLocBas];
    double * volctrl_z = new double [nLocBas];

    int * rotated_local_ien = new int [nLocBas];
    double * rotated_local_disp = new double [nLocBas * 3];

    int rotated_face_id = -1;
    int rotated_tag = tag;
    int num_rotated_ele = itf_part->get_num_tagged_rotated_ele(itf_id, rotated_tag);
    // SYS_T::commPrint("    num_rotated_ele:%d\n", num_rotated_ele);

    for(int ee_index{0}; ee_index<num_rotated_ele; ++ee_index)
    {
      const int ee = itf_part->get_tagged_rotated_ele(itf_id, rotated_tag, ee_index);
      itf_part->get_rotated_ele_ctrlPts(itf_id, ee, inictrl_x, inictrl_y, inictrl_z);
      SI_sol->get_rotated_mdisp(itf_part, itf_id, ee, rotated_local_ien, rotated_local_disp);
      SI_T::get_currPts(inictrl_x, inictrl_y, inictrl_z, rotated_local_disp, nLocBas, volctrl_x, volctrl_y, volctrl_z);
      
      rotated_face_id = itf_part->get_rotated_face_id(itf_id, ee);

      const auto facectrl = rotated_elementv->get_face_ctrlPts(rotated_face_id,
        volctrl_x, volctrl_y, volctrl_z);

      rotated_xi->reset();
      is_found = FE_T::search_closest_point(fixed_pt, elements,
        facectrl[0].data(), facectrl[1].data(), facectrl[2].data(), rotated_xi);

      if(is_found)
      {
        rotated_ee = ee;
        // SYS_T::commPrint("  found in rotated_ee = %d.\n\n", rotated_ee);
        break;
      }
    }

    // Second try
    if(is_found == false && tag != 0)
    {
      rotated_tag = tag - 1;

      num_rotated_ele = itf_part->get_num_tagged_rotated_ele(itf_id, rotated_tag);
      // SYS_T::commPrint("    num_rotated_ele:%d\n", num_rotated_ele);

      for(int ee_index{0}; ee_index<num_rotated_ele; ++ee_index)
      {
        const int ee = itf_part->get_tagged_rotated_ele(itf_id, rotated_tag, ee_index);
        itf_part->get_rotated_ele_ctrlPts(itf_id, ee, inictrl_x, inictrl_y, inictrl_z);
        SI_sol->get_rotated_mdisp(itf_part, itf_id, ee, rotated_local_ien, rotated_local_disp);
        SI_T::get_currPts(inictrl_x, inictrl_y, inictrl_z, rotated_local_disp, nLocBas, volctrl_x, volctrl_y, volctrl_z);
        
        rotated_face_id = itf_part->get_rotated_face_id(itf_id, ee);

        const auto facectrl = rotated_elementv->get_face_ctrlPts(rotated_face_id,
          volctrl_x, volctrl_y, volctrl_z);

        rotated_xi->reset();
        is_found = FE_T::search_closest_point(fixed_pt, elements,
          facectrl[0].data(), facectrl[1].data(), facectrl[2].data(), rotated_xi);

        if(is_found)
        {
          rotated_ee = ee;
          // SYS_T::commPrint("  found in rotated_ee = %d.\n\n", rotated_ee);
          break;
        }
      }
    }

    // Third try
    if(is_found == false && tag != itf_part->get_num_tag(itf_id) - 1)
    {
      rotated_tag = tag + 1;

      num_rotated_ele = itf_part->get_num_tagged_rotated_ele(itf_id, rotated_tag);
      // SYS_T::commPrint("    num_rotated_ele:%d\n", num_rotated_ele);

      for(int ee_index{0}; ee_index<num_rotated_ele; ++ee_index)
      {
        const int ee = itf_part->get_tagged_rotated_ele(itf_id, rotated_tag, ee_index);
        itf_part->get_rotated_ele_ctrlPts(itf_id, ee, inictrl_x, inictrl_y, inictrl_z);
        SI_sol->get_rotated_mdisp(itf_part, itf_id, ee, rotated_local_ien, rotated_local_disp);
        SI_T::get_currPts(inictrl_x, inictrl_y, inictrl_z, rotated_local_disp, nLocBas, volctrl_x, volctrl_y, volctrl_z);
        
        rotated_face_id = itf_part->get_rotated_face_id(itf_id, ee);

        const auto facectrl = rotated_elementv->get_face_ctrlPts(rotated_face_id,
          volctrl_x, volctrl_y, volctrl_z);

        rotated_xi->reset();
        is_found = FE_T::search_closest_point(fixed_pt, elements,
          facectrl[0].data(), facectrl[1].data(), facectrl[2].data(), rotated_xi);

        if(is_found)
        {
          rotated_ee = ee;
          // SYS_T::commPrint("  found in rotated_ee = %d.\n\n", rotated_ee);
          break;
        }
      }
    }

    delete [] inictrl_x; inictrl_x = nullptr;
    delete [] inictrl_y; inictrl_y = nullptr;
    delete [] inictrl_z; inictrl_z = nullptr;

    delete [] volctrl_x; volctrl_x = nullptr;
    delete [] volctrl_y; volctrl_y = nullptr;
    delete [] volctrl_z; volctrl_z = nullptr;

    delete [] rotated_local_ien; rotated_local_ien = nullptr;
    delete [] rotated_local_disp; rotated_local_disp = nullptr;

    SYS_T::print_fatal_if(is_found == false,
      "Error, SI_quad_point::search_opposite_point: cannot find opposite rotated point.\n");

    // tag = rotated_tag;
  }

  void SI_quad_point::search_opposite_fixed_point(
    const Vector_3 &rotated_pt,
    const ALocal_Interface * const &itf_part,
    const SI_solution * const &SI_sol,
    const int &itf_id,
    FEAElement * fixed_elementv,
    FEAElement * elements,
    int &tag,
    int &fixed_ee,
    IQuadPts * const &fixed_xi )
  {
    const int nLocBas = itf_part->get_nLocBas();
    bool is_found = false;

    double * volctrl_x = new double [nLocBas];
    double * volctrl_y = new double [nLocBas];
    double * volctrl_z = new double [nLocBas];

    int * fixed_local_ien = new int [nLocBas];

    int fixed_face_id = -1;
    int fixed_tag = tag;
    int num_fixed_ele = itf_part->get_num_tagged_fixed_ele(itf_id, fixed_tag);

    for(int ee_index{0}; ee_index<num_fixed_ele; ++ee_index)
    {
      const int ee = itf_part->get_tagged_fixed_ele(itf_id, fixed_tag, ee_index);
      itf_part->get_fixed_ele_ctrlPts(itf_id, ee, volctrl_x, volctrl_y, volctrl_z);
      
      fixed_face_id = itf_part->get_fixed_face_id(itf_id, ee);

      const auto facectrl = fixed_elementv->get_face_ctrlPts(fixed_face_id,
        volctrl_x, volctrl_y, volctrl_z);

      fixed_xi->reset();
      is_found = FE_T::search_closest_point(rotated_pt, elements,
        facectrl[0].data(), facectrl[1].data(), facectrl[2].data(), fixed_xi);

      if(is_found)
      {
        fixed_ee = ee;
        break;
      }
    }

    // Second try
    if(is_found == false && tag != 0)
    {
      fixed_tag = tag - 1;

      num_fixed_ele = itf_part->get_num_tagged_fixed_ele(itf_id, fixed_tag);

      for(int ee_index{0}; ee_index<num_fixed_ele; ++ee_index)
      {
        const int ee = itf_part->get_tagged_fixed_ele(itf_id, fixed_tag, ee_index);
        itf_part->get_fixed_ele_ctrlPts(itf_id, ee, volctrl_x, volctrl_y, volctrl_z);
        
        fixed_face_id = itf_part->get_fixed_face_id(itf_id, ee);

        const auto facectrl = fixed_elementv->get_face_ctrlPts(fixed_face_id,
          volctrl_x, volctrl_y, volctrl_z);

        fixed_xi->reset();
        is_found = FE_T::search_closest_point(rotated_pt, elements,
          facectrl[0].data(), facectrl[1].data(), facectrl[2].data(), fixed_xi);

        if(is_found)
        {
          fixed_ee = ee;
          break;
        }
      }
    }

    // Third try
    if(is_found == false && tag != itf_part->get_num_tag(itf_id) - 1)
    {
      fixed_tag = tag + 1;

      num_fixed_ele = itf_part->get_num_tagged_fixed_ele(itf_id, fixed_tag);

      for(int ee_index{0}; ee_index<num_fixed_ele; ++ee_index)
      {
        const int ee = itf_part->get_tagged_fixed_ele(itf_id, fixed_tag, ee_index);
        itf_part->get_fixed_ele_ctrlPts(itf_id, ee, volctrl_x, volctrl_y, volctrl_z);
        
        fixed_face_id = itf_part->get_fixed_face_id(itf_id, ee);

        const auto facectrl = fixed_elementv->get_face_ctrlPts(fixed_face_id,
          volctrl_x, volctrl_y, volctrl_z);

        fixed_xi->reset();
        is_found = FE_T::search_closest_point(rotated_pt, elements,
          facectrl[0].data(), facectrl[1].data(), facectrl[2].data(), fixed_xi);

        if(is_found)
        {
          fixed_ee = ee;
          break;
        }
      }
    }

    delete [] volctrl_x; volctrl_x = nullptr;
    delete [] volctrl_y; volctrl_y = nullptr;
    delete [] volctrl_z; volctrl_z = nullptr;

    delete [] fixed_local_ien; fixed_local_ien = nullptr;

    SYS_T::print_fatal_if(is_found == false,
      "Error, SI_quad_point::search_opposite_point: cannot find opposite fixed point.\n");
  }

  void get_currPts( const double * const &ept_x,
    const double * const &ept_y,
    const double * const &ept_z,
    const double * const &disp,
    const int &len,
    double * const &currPt_x,
    double * const &currPt_y,
    double * const &currPt_z )
    {
      for(int ii=0; ii<len; ++ii)
      {
        currPt_x[ii] = ept_x[ii] + disp[3*ii];
        currPt_y[ii] = ept_y[ii] + disp[3*ii+1];
        currPt_z[ii] = ept_z[ii] + disp[3*ii+2];
      }
    }
}

// EOF
