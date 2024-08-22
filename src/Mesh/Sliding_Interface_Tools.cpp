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
    fixed_node_part_tag.resize(num_itf);
    fixed_node_loc_pos.resize(num_itf);

    num_rotated_node.assign(num_itf, 0);
    rotated_node_sol.resize(num_itf);
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

      fixed_node_part_tag[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_node_part_tag" );

      fixed_node_loc_pos[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_node_loc_pos" );

      num_rotated_node[ii] = VEC_T::get_size(h5r -> read_intVector( subgroup_name.c_str(), "rotated_node_map" ));

      rotated_node_sol[ii] = std::vector<double> (dof_sol * num_rotated_node[ii], 0.0);

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


  SI_quad_point::SI_quad_point(const ALocal_Interface * const &itf, const int &nqp_sur_in)
  {
    nqp_sur = nqp_sur_in;
    curr_tag.resize(itf->get_num_itf());
    curr_ee.resize(itf->get_num_itf());
    curr_xi.resize(itf->get_num_itf());

    for(int ii = 0; ii < itf->get_num_itf(); ++ii)
    {
      curr_tag[ii].assign(itf->get_num_fixed_ele(ii) * nqp_sur, -1);
      curr_ee[ii].assign(itf->get_num_fixed_ele(ii) * nqp_sur, -1);
      curr_xi[ii].assign(itf->get_num_fixed_ele(ii) * nqp_sur, std::vector<double>(2, 0.0));
    }
  }

  void SI_quad_point::set_curr(const double &itf_id, const int &fixed_ee, const int &qua,
    const int &ele_tag, const int &rotated_ee, const std::vector<double> &xi)
  {
    curr_tag[itf_id][fixed_ee * nqp_sur + qua] = ele_tag;
    curr_ee[itf_id][fixed_ee * nqp_sur + qua] = rotated_ee;
    curr_xi[itf_id][fixed_ee * nqp_sur + qua][0] = xi[0];
    curr_xi[itf_id][fixed_ee * nqp_sur + qua][1] = xi[1];
  }

  void SI_quad_point::get_curr(const double &itf_id, const int &fixed_ee, const int &qua,
    int &ele_tag, int &rotated_ee, std::vector<double> &xi) const
  {
    ele_tag = curr_tag[itf_id][fixed_ee * nqp_sur + qua];
    rotated_ee = curr_ee[itf_id][fixed_ee * nqp_sur + qua];
    xi = curr_xi[itf_id][fixed_ee * nqp_sur + qua];
  }

  void SI_quad_point::search_all_opposite_point(
    const double &curr_time,
    FEAElement * const &fixed_elementv,
    FEAElement * const &rotated_elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_s,
    IQuadPts * const &free_quad,
    ALocal_Interface * const &itf_part )
  {
    const int nLocBas = itf_part->get_nLocBas();
    double * ctrl_x = new double [nLocBas];
    double * ctrl_y = new double [nLocBas];
    double * ctrl_z = new double [nLocBas];

    int * fixed_local_ien = new int [nLocBas];

    int * rotated_local_ien = new int [nLocBas];

    const int num_itf {itf_part->get_num_itf()};

    for(int itf_id{0}; itf_id<num_itf; ++itf_id)
    {
      SYS_T::commPrint("itf_id = %d\n", itf_id);
      const int num_fixed_elem = itf_part->get_num_fixed_ele(itf_id);

      for(int ee{0}; ee<num_fixed_elem; ++ee)
      {
        // SYS_T::commPrint("  fixed_ee = %d\n", ee);
        // const int local_ee_index{itf_part->get_fixed_ele_id(itf_id, ee)};

        itf_part->get_fixed_ele_ctrlPts(itf_id, ee, ctrl_x, ctrl_y, ctrl_z);

        const int fixed_face_id {itf_part->get_fixed_face_id(itf_id, ee)};

        int ele_tag {itf_part->get_fixed_ele_tag(itf_id, ee)};

        fixed_elementv->buildBasis(fixed_face_id, quad_s, ctrl_x, ctrl_y, ctrl_z);

        const int fixed_face_nqp {quad_s->get_num_quadPts()};

        std::vector<double> R(nLocBas, 0.0);

        for(int qua{0}; qua<fixed_face_nqp; ++qua)
        {
          fixed_elementv->get_R(qua, &R[0]);

          // The xyz-coordinates of the quadrature point
          Vector_3 coor(0.0, 0.0, 0.0);
          for(int ii{0}; ii<nLocBas; ++ii)
          {
            coor.x() += ctrl_x[ii] * R[ii];
            coor.y() += ctrl_y[ii] * R[ii];
            coor.z() += ctrl_z[ii] * R[ii];
          }

          // SYS_T::commPrint("    point %d:\n", qua);

          int rotated_ee {0};
          search_opposite_point(curr_time, coor, itf_part, itf_id, rotated_elementv, elements, ele_tag, rotated_ee, free_quad);

          std::vector<double> rotated_xi = {free_quad->get_qp(0, 0), free_quad->get_qp(0, 1)};

          set_curr(itf_id, ee, qua, ele_tag, rotated_ee, rotated_xi);
        }
      }
    }

    delete [] fixed_local_ien; fixed_local_ien = nullptr;
    delete [] rotated_local_ien; rotated_local_ien = nullptr;

    delete [] ctrl_x; ctrl_x = nullptr;
    delete [] ctrl_y; ctrl_y = nullptr;
    delete [] ctrl_z; ctrl_z = nullptr;
  }
  
  void SI_quad_point::search_opposite_point(
    const double &curr_time,
    const Vector_3 &fixed_pt,
    const ALocal_Interface * const &itf_part,
    const int &itf_id,
    FEAElement * rotated_elementv,
    FEAElement * elements,
    int &tag,
    int &rotated_ee,
    IQuadPts * const &rotated_xi )
  {
    const int nLocBas = itf_part->get_nLocBas();
    bool is_found = false;

    double * volctrl_x = new double [nLocBas];
    double * volctrl_y = new double [nLocBas];
    double * volctrl_z = new double [nLocBas];

    const int snlocbas = elements->get_nLocBas();

    int rotated_face_id = -1;
    int rotated_tag = tag;
    int num_rotated_ele = itf_part->get_num_rotated_ele(itf_id, rotated_tag);
    // SYS_T::commPrint("    num_rotated_ele:%d\n", num_rotated_ele);

    for(int ee{0}; ee<num_rotated_ele; ++ee)
    {
      // SYS_T::commPrint("    search rotated ee = %d\n", ee);
      itf_part->get_rotated_ele_ctrlPts(itf_id, rotated_tag, ee, curr_time, volctrl_x, volctrl_y, volctrl_z);
      
      rotated_face_id = itf_part->get_rotated_face_id(itf_id, rotated_tag, ee);

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

      num_rotated_ele = itf_part->get_num_rotated_ele(itf_id, rotated_tag);
      // SYS_T::commPrint("    num_rotated_ele:%d\n", num_rotated_ele);

      for(int ee{0}; ee<num_rotated_ele; ++ee)
      {
        // SYS_T::commPrint("    search rotated ee = %d\n", ee);
        itf_part->get_rotated_ele_ctrlPts(itf_id, rotated_tag, ee, curr_time, volctrl_x, volctrl_y, volctrl_z);
        
        rotated_face_id = itf_part->get_rotated_face_id(itf_id, rotated_tag, ee);

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

      num_rotated_ele = itf_part->get_num_rotated_ele(itf_id, rotated_tag);
      // SYS_T::commPrint("    num_rotated_ele:%d\n", num_rotated_ele);

      for(int ee{0}; ee<num_rotated_ele; ++ee)
      {
        // SYS_T::commPrint("    search rotated ee = %d\n", ee);
        itf_part->get_rotated_ele_ctrlPts(itf_id, rotated_tag, ee, curr_time, volctrl_x, volctrl_y, volctrl_z);
        
        rotated_face_id = itf_part->get_rotated_face_id(itf_id, rotated_tag, ee);

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

    delete [] volctrl_x; volctrl_x = nullptr;
    delete [] volctrl_y; volctrl_y = nullptr;
    delete [] volctrl_z; volctrl_z = nullptr;

    SYS_T::print_fatal_if(is_found == false,
      "Error, SI_quad_point::search_opposite_point: cannot find opposite point.\n");

    tag = rotated_tag;
  }
}

// EOF
