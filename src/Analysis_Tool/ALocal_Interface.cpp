#include "ALocal_Interface.hpp"

ALocal_Interface::ALocal_Interface( const std::string &fileBaseName, const int &cpu_rank)
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/sliding");

  num_itf = h5r -> read_intScalar( gname.c_str(), "num_interface" );
  SYS_T::print_fatal_if(num_itf < 1, "Error, ALocal_Interface: there is no interface in this geometric model.\n");

  num_fixed_ele = h5r -> read_intVector( gname.c_str(), "num_part_fixed_cell" );
  num_tag = h5r -> read_intVector( gname.c_str(), "num_tag" );

  std::string groupbase(gname);
  groupbase.append("/interfaceid_");

  num_fixed_node.assign(num_itf, 0);
  num_rotated_node.assign(num_itf, 0);

  num_rotated_ele.resize(num_itf);

  fixed_ele_face_id.resize(num_itf);
  fixed_ele_tag.resize(num_itf);
  fixed_layer_ien.resize(num_itf);
  fixed_node_xyz.resize(num_itf);
  fixed_node_sol.resize(num_itf);
  fixed_node_part_tag.resize(num_itf);
  fixed_node_loc_pos.resize(num_itf);
  fixed_node_id.resize(num_itf);
  fixed_ID.resize(num_itf);

  rotated_layer_ien.resize(num_itf);
  rotated_layer_face_id.resize(num_itf);
  init_rotated_node_xyz.resize(num_itf);
  rotated_node_sol.resize(num_itf);
  rotated_node_mvelo.resize(num_itf);
  rotated_node_disp.resize(num_itf); 
  rotated_node_part_tag.resize(num_itf);
  rotated_node_loc_pos.resize(num_itf);
  rotated_node_id.resize(num_itf);
  rotated_ID.resize(num_itf);

  const std::string mesh_info("/Global_Mesh_Info");
  nLocBas = h5r -> read_intScalar(mesh_info.c_str(), "nLocBas");
  dof_sol = h5r -> read_intScalar(mesh_info.c_str(), "dofNum");

  const std::string part_node("/Local_Node");
  nlgn = h5r -> read_intScalar(part_node.c_str(), "nlocghonode");

  cpu = cpu_rank;

  for(int ii=0; ii<num_itf; ++ii)
  {
    std::string subgroup_name(groupbase);
    subgroup_name.append( std::to_string(ii) );

    fixed_ele_face_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_cell_face_id" );

    fixed_ele_tag[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_cell_tag" );

    fixed_layer_ien[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_cell_ien" );

    fixed_node_xyz[ii] = h5r -> read_doubleVector( subgroup_name.c_str(), "fixed_node_xyz" );

    fixed_node_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_node_map" );

    fixed_ID[ii] = h5r -> read_intVector(  subgroup_name.c_str(), "fixed_ID" );

    num_fixed_node[ii] = VEC_T::get_size(fixed_node_id[ii]);

    fixed_node_sol[ii] = std::vector<double> (dof_sol * num_fixed_node[ii], 0.0);

    fixed_node_part_tag[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_node_part_tag" );

    fixed_node_loc_pos[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_node_loc_pos" );

    rotated_layer_ien[ii].resize(num_tag[ii]);
    rotated_layer_face_id[ii].resize(num_tag[ii]);
    num_rotated_ele[ii].resize(num_tag[ii]);
    
    std::string subsubgroupbase(subgroup_name);
    subsubgroupbase.append("/tag_");

    for(int jj=0; jj<num_tag[ii]; ++jj)
    {
      std::string subsubgroup_name(subsubgroupbase);
      subsubgroup_name.append( std::to_string(jj) );

      num_rotated_ele[ii][jj] = h5r -> read_intScalar( subsubgroup_name.c_str(), "num_rotated_cell" );

      if(num_rotated_ele[ii][jj] > 0)
      {
        rotated_layer_ien[ii][jj] = h5r -> read_intVector( subsubgroup_name.c_str(), "rotated_cell_ien" );

        rotated_layer_face_id[ii][jj] = h5r -> read_intVector( subsubgroup_name.c_str(), "rotated_cell_face_id" );
      }
    }

    init_rotated_node_xyz[ii] = h5r -> read_doubleVector( subgroup_name.c_str(), "rotated_node_xyz" );

    rotated_node_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "rotated_node_map" );

    rotated_ID[ii] = h5r -> read_intVector( subgroup_name.c_str(), "rotated_ID" );

    num_rotated_node[ii] = VEC_T::get_size(rotated_node_id[ii]);

    rotated_node_sol[ii] = std::vector<double> (dof_sol * num_rotated_node[ii], 0.0);

    rotated_node_mvelo[ii] = std::vector<double> (dof_sol * num_rotated_node[ii], 0.0);

    rotated_node_disp[ii] = std::vector<double> (dof_sol * num_rotated_node[ii], 0.0); 

    rotated_node_part_tag[ii] = h5r -> read_intVector( subgroup_name.c_str(), "rotated_node_part_tag" );

    rotated_node_loc_pos[ii] = h5r -> read_intVector( subgroup_name.c_str(), "rotated_node_loc_pos" );
  }

  delete h5r; H5Fclose( file_id );
}

void ALocal_Interface::print_info() const
{
  SYS_T::commPrint("Interfaces: %d\n", num_itf);
}

void ALocal_Interface::get_fixed_ele_ctrlPts(const int &ii, const int &ee,
  double * const volctrl_x,  double * const volctrl_y,  double * const volctrl_z) const
{
  for(int nn{0}; nn < nLocBas; ++nn)
  {
    int node = get_fixed_layer_ien(ii, nLocBas * ee + nn);

    volctrl_x[nn] = get_fixed_node_xyz(ii, 3 * node);
    volctrl_y[nn] = get_fixed_node_xyz(ii, 3 * node + 1);
    volctrl_z[nn] = get_fixed_node_xyz(ii, 3 * node + 2);
  }
}

void ALocal_Interface::get_rotated_ele_ctrlPts(const int &ii, const int &tag, const int &ee, const double &tt,
  double * const volctrl_x, double * const volctrl_y, double * const volctrl_z) const
{
  for(int nn{0}; nn < nLocBas; ++nn)
  {
    int node = get_rotated_layer_ien(ii, tag, nLocBas * ee + nn);

    volctrl_x[nn] = get_init_rotated_node_xyz(ii, 3 * node);
    volctrl_y[nn] = get_init_rotated_node_xyz(ii, 3 * node + 1);
    volctrl_z[nn] = get_init_rotated_node_xyz(ii, 3 * node + 2);
  }
}

void ALocal_Interface::restore_node_sol(const PDNSolution * const &sol)
{
  double * array = new double [nlgn * dof_sol];

  sol->GetLocalArray( array );

  Zero_node_sol();

  for(int ii = 0; ii < num_itf; ++ii)
  { 
    std::vector<double> temp_fixed_node_sol = fixed_node_sol[ii];
    for(int nn = 0; nn < num_fixed_node[ii]; ++nn)
    {
      // Pick out the solution of nodes in each part
      if(fixed_node_part_tag[ii][nn] == cpu)
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
      if(rotated_node_part_tag[ii][nn] == cpu)
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

  delete [] array;
}

void ALocal_Interface::restore_node_mvelo(const PDNSolution * const &mvelo)
{
  double * array = new double [nlgn * 3];

  mvelo->GetLocalArray( array );

  Zero_node_mvelo();

  for(int ii = 0; ii < num_itf; ++ii)
  { 
    std::vector<double> temp_rotated_node_mvelo = rotated_node_mvelo[ii];
    for(int nn = 0; nn < num_rotated_node[ii]; ++nn)
    {
      // Pick out the solution of nodes in each part
      if(rotated_node_part_tag[ii][nn] == cpu)
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

  delete [] array;  
}

void ALocal_Interface::restore_node_disp(const PDNSolution * const &disp)
{
  double * array = new double [nlgn * 3];

  disp->GetLocalArray( array );

  Zero_node_disp();

  for(int ii = 0; ii < num_itf; ++ii)
  { 
    std::vector<double> temp_rotated_node_disp = rotated_node_disp[ii];
    for(int nn = 0; nn < num_rotated_node[ii]; ++nn)
    {
      // Pick out the solution of nodes in each part
      if(rotated_node_part_tag[ii][nn] == cpu)
      { 
        // Like GetLocal in PGAssem
        const int loc_pos = rotated_node_loc_pos[ii][nn];
        for(int dd = 0; dd < 3; ++dd)
          temp_rotated_node_disp[3 * nn + dd] = array[3 * loc_pos + dd];
      }
    }

    // Summation from each part
    MPI_Allreduce(&temp_rotated_node_disp[0], &rotated_node_disp[ii][0], 3 * num_rotated_node[ii], MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  }

  delete [] array;  
}

void ALocal_Interface::init_curr(const int &nqp_sur_in)
{
  nqp_sur = nqp_sur_in;
  curr_tag.resize(num_itf);
  curr_ee.resize(num_itf);
  curr_xi.resize(num_itf);

  for(int ii = 0; ii < num_itf; ++ii)
  {
    curr_tag[ii].assign(num_fixed_ele[ii] * nqp_sur, -1);
    curr_ee[ii].assign(num_fixed_ele[ii] * nqp_sur, -1);
    curr_xi[ii].assign(num_fixed_ele[ii] * nqp_sur, std::vector<double>(2, 0.0));
  }
}

void ALocal_Interface::set_curr(const double &itf_id, const int &fixed_ee, const int &qua,
  const int &ele_tag, const int &rotated_ee, const std::vector<double> &xi)
{
  curr_tag[itf_id][fixed_ee * nqp_sur + qua] = ele_tag;
  curr_ee[itf_id][fixed_ee * nqp_sur + qua] = rotated_ee;
  curr_xi[itf_id][fixed_ee * nqp_sur + qua][0] = xi[0];
  curr_xi[itf_id][fixed_ee * nqp_sur + qua][1] = xi[1];
}

void ALocal_Interface::get_curr(const double &itf_id, const int &fixed_ee, const int &qua,
  int &ele_tag, int &rotated_ee, std::vector<double> &xi) const
{
  ele_tag = curr_tag[itf_id][fixed_ee * nqp_sur + qua];
  rotated_ee = curr_ee[itf_id][fixed_ee * nqp_sur + qua];
  xi = curr_xi[itf_id][fixed_ee * nqp_sur + qua];
}

// EOF
