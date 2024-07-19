#include "Interface_pair.hpp"

Interface_pair::Interface_pair(const std::string &fixed_vtkfile, const std::string &rotated_vtkfile,
  const std::string &fixed_h5file, const int &total_num_fixed_elem, const int &total_num_fixed_pt,
  const std::vector<double> &all_vol_ctrlPts, const IIEN * const &VIEN, const int &elemtype_in,
  const std::vector<double> &intervals_in, const int &direction_in) :
  interface_type {0}, T0_axial_direction{direction_in}, T1_surface_centroid{Vector_3(0,0,0)}
{
  Initialize(fixed_vtkfile, rotated_vtkfile, fixed_h5file, total_num_fixed_elem, total_num_fixed_pt,
    all_vol_ctrlPts, VIEN, elemtype_in, intervals_in);
}

Interface_pair::Interface_pair(const std::string &fixed_vtkfile, const std::string &rotated_vtkfile,
  const std::string &fixed_h5file, const int &total_num_fixed_elem, const int &total_num_fixed_pt,
  const std::vector<double> &all_vol_ctrlPts, const IIEN * const &VIEN, const int &elemtype_in,
  const std::vector<double> &intervals_in, const Vector_3 &centroid_in) :
  interface_type {1}, T0_axial_direction{-1}, T1_surface_centroid{centroid_in}
{
  Initialize(fixed_vtkfile, rotated_vtkfile, fixed_h5file, total_num_fixed_elem, total_num_fixed_pt,
    all_vol_ctrlPts, VIEN, elemtype_in, intervals_in);
}

void Interface_pair::Initialize(const std::string &fixed_vtkfile,
                    const std::string &rotated_vtkfile,
                    const std::string &fixed_h5file,
                    const int &total_num_fixed_elem,
                    const int &total_num_fixed_pt,
                    const std::vector<double> &all_vol_ctrlPts,
                    const IIEN * const &VIEN,
                    const int &elemtype_in,
                    const std::vector<double> &intervals_in)
{ 
  // Read the vtk files
  int num_fixed_sur_node {0};
  std::vector<double> fixed_sur_pt_xyz {};
  std::vector<int> fixed_sur_ien {};

  VTK_T::read_grid(fixed_vtkfile, num_fixed_sur_node, num_fixed_ele, fixed_sur_pt_xyz, fixed_sur_ien);

  std::vector<int> fixed_global_node = VTK_T::read_int_PointData(fixed_vtkfile, "GlobalNodeID");
  std::vector<int> fixed_global_cell = VTK_T::read_int_CellData(fixed_vtkfile, "GlobalElementID");

  int num_rotated_sur_node {0};
  std::vector<double> rotated_sur_pt_xyz {};
  std::vector<int> rotated_sur_ien {};
  VTK_T::read_grid(rotated_vtkfile, num_rotated_sur_node, num_rotated_ele, rotated_sur_pt_xyz, rotated_sur_ien);

  std::vector<int> rotated_global_node = VTK_T::read_int_PointData(rotated_vtkfile, "GlobalNodeID");
  std::vector<int> rotated_global_cell = VTK_T::read_int_CellData(rotated_vtkfile, "GlobalElementID");

  switch (elemtype_in)
  {
    case 501:
    {
      s_nLocBas = 3; v_nLocBas = 4;
    }
    break;

    case 502:
    {
      s_nLocBas = 6; v_nLocBas = 10;
    }
    break;

    case 601:
    {
      s_nLocBas = 6; v_nLocBas = 8;
    }
    break;

    case 602:
    {
      s_nLocBas = 9; v_nLocBas = 27;
    }
    break;
    
    default:
      SYS_T::print_fatal("Error, Interface_pair: unknown element type.\n");
    break;
  }

  fixed_face_id.resize(num_fixed_ele);
  fixed_vien.resize(v_nLocBas * num_fixed_ele);

  rotated_face_id.resize(num_rotated_ele);
  rotated_vien.resize(v_nLocBas * num_rotated_ele);
    
  // Read the partion tag from the .h5 file
  fixed_cpu_rank = HDF5_T::read_intVector( fixed_h5file.c_str(), "/", "part");

  // Generate the face id and layer's ien array
  if(elemtype_in == 501 || elemtype_in == 502)
  {
    TET_T::Tet4 * tetcell = new TET_T::Tet4();
    
    for(int ee=0; ee<num_fixed_ele; ++ee)
    {
      const int node_t[3] {fixed_sur_ien[ee * s_nLocBas + 0],
                            fixed_sur_ien[ee * s_nLocBas + 1],
                            fixed_sur_ien[ee * s_nLocBas + 2]};

      const int node_t_gi[3] {fixed_global_node[node_t[0]],
                              fixed_global_node[node_t[1]],
                              fixed_global_node[node_t[2]]};
      
      const int cell_gi = fixed_global_cell[ee];

      const int tet_n[4] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                            VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3) };

      tetcell->reset(tet_n[0], tet_n[1], tet_n[2], tet_n[3]);
      fixed_face_id[ee] = tetcell->get_face_id(node_t_gi[0], node_t_gi[1], node_t_gi[2]);

      for(int ii=0; ii<v_nLocBas; ++ii)
        fixed_vien[ee * v_nLocBas + ii] = VIEN->get_IEN(cell_gi, ii);
    }

    for(int ee=0; ee<num_rotated_ele; ++ee)
    {
      const int node_t[3] {rotated_sur_ien[ee * s_nLocBas + 0],
                            rotated_sur_ien[ee * s_nLocBas + 1],
                            rotated_sur_ien[ee * s_nLocBas + 2]};

      const int node_t_gi[3] {rotated_global_node[node_t[0]] + total_num_fixed_pt,
                              rotated_global_node[node_t[1]] + total_num_fixed_pt,
                              rotated_global_node[node_t[2]] + total_num_fixed_pt};

      const int cell_gi = rotated_global_cell[ee] + total_num_fixed_elem;

      const int tet_n[4] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                            VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3) };

      tetcell->reset(tet_n[0], tet_n[1], tet_n[2], tet_n[3]);
      rotated_face_id[ee] = tetcell->get_face_id(node_t_gi[0], node_t_gi[1], node_t_gi[2]);

      for(int ii=0; ii<v_nLocBas; ++ii)
        rotated_vien[ee * v_nLocBas + ii] = VIEN->get_IEN(cell_gi, ii);
    }

    delete tetcell;
  }
  else if(elemtype_in == 601 || elemtype_in == 602)
  {
    HEX_T::Hex8 * hexcell = new HEX_T::Hex8();

    for(int ee=0; ee<num_fixed_ele; ++ee)
    {
      const int node_q[4] {fixed_sur_ien[ee * s_nLocBas + 0],
                           fixed_sur_ien[ee * s_nLocBas + 1],
                           fixed_sur_ien[ee * s_nLocBas + 2],
                           fixed_sur_ien[ee * s_nLocBas + 3]};

      const int node_q_gi[4] {fixed_global_node[node_q[0]],
                              fixed_global_node[node_q[1]],
                              fixed_global_node[node_q[2]],
                              fixed_global_node[node_q[3]]};
      
      const int cell_gi = fixed_global_cell[ee];

      const int hex_n[8] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                           VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3),
                           VIEN->get_IEN(cell_gi, 5), VIEN->get_IEN(cell_gi, 6),
                           VIEN->get_IEN(cell_gi, 7), VIEN->get_IEN(cell_gi, 8), };

      hexcell->reset(hex_n[0], hex_n[1], hex_n[2], hex_n[3],
                      hex_n[4], hex_n[5], hex_n[6], hex_n[7]);
      fixed_face_id[ee] = hexcell->get_face_id(node_q_gi[0], node_q_gi[1], node_q_gi[2], node_q_gi[3]);

      for(int ii=0; ii<v_nLocBas; ++ii)
        fixed_vien[ee * v_nLocBas + ii] = VIEN->get_IEN(cell_gi, ii);
    }

    for(int ee=0; ee<num_rotated_ele; ++ee)
    {
      const int node_q[4] {rotated_sur_ien[ee * s_nLocBas + 0],
                            rotated_sur_ien[ee * s_nLocBas + 1],
                            rotated_sur_ien[ee * s_nLocBas + 2],
                            rotated_sur_ien[ee * s_nLocBas + 3]};

      const int node_q_gi[4] {rotated_global_node[node_q[0]] + total_num_fixed_pt,
                              rotated_global_node[node_q[1]] + total_num_fixed_pt,
                              rotated_global_node[node_q[2]] + total_num_fixed_pt,
                              rotated_global_node[node_q[3]] + total_num_fixed_pt};

      const int cell_gi = rotated_global_cell[ee] + total_num_fixed_elem;

      const int hex_n[8] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                            VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3),
                            VIEN->get_IEN(cell_gi, 5), VIEN->get_IEN(cell_gi, 6),
                            VIEN->get_IEN(cell_gi, 7), VIEN->get_IEN(cell_gi, 8), };

      hexcell->reset(hex_n[0], hex_n[1], hex_n[2], hex_n[3],
                      hex_n[4], hex_n[5], hex_n[6], hex_n[7]);
      rotated_face_id[ee] = hexcell->get_face_id(node_q_gi[0], node_q_gi[1], node_q_gi[2], node_q_gi[3]);

      for(int ii=0; ii<v_nLocBas; ++ii)
        rotated_vien[ee * v_nLocBas + ii] = VIEN->get_IEN(cell_gi, ii);
    }

    delete hexcell;
  }
  else
    SYS_T::print_fatal("Error: Interface_pair, unknown element type.\n");

  // Generate the global node id and xyz
  fixed_global_node = fixed_vien;
  VEC_T::sort_unique_resize(fixed_global_node);
  const int num_fixed_node = VEC_T::get_size(fixed_global_node);

  // PERIGEE_OMP_PARALLEL_FOR
  for(int &nodeid : fixed_vien)
  {
    const int local_id = VEC_T::get_pos(fixed_global_node, nodeid);
    nodeid = local_id;
  }

  fixed_pt_xyz.resize(3 * num_fixed_node);
  
  PERIGEE_OMP_PARALLEL_FOR
  for(int nn=0; nn<num_fixed_node; ++nn)
  {
    const int GID = fixed_global_node[nn];
    fixed_pt_xyz[3 * nn]     = all_vol_ctrlPts[3 * GID];
    fixed_pt_xyz[3 * nn + 1] = all_vol_ctrlPts[3 * GID + 1];
    fixed_pt_xyz[3 * nn + 2] = all_vol_ctrlPts[3 * GID + 2];
  }

  rotated_global_node = rotated_vien;
  VEC_T::sort_unique_resize(rotated_global_node);
  const int num_rotated_node = VEC_T::get_size(rotated_global_node);

  // PERIGEE_OMP_PARALLEL_FOR
  for(int &nodeid : rotated_vien)
  {
    const int local_id = VEC_T::get_pos(rotated_global_node, nodeid);
    nodeid = local_id;
  }

  rotated_pt_xyz.resize(3 * num_rotated_node);
  
  PERIGEE_OMP_PARALLEL_FOR
  for(int nn=0; nn<num_rotated_node; ++nn)
  {
    const int GID = rotated_global_node[nn];
    rotated_pt_xyz[3 * nn]     = all_vol_ctrlPts[3 * GID];
    rotated_pt_xyz[3 * nn + 1] = all_vol_ctrlPts[3 * GID + 1];
    rotated_pt_xyz[3 * nn + 2] = all_vol_ctrlPts[3 * GID + 2];
  }

  Check_interval(intervals_in);

  // Generate the interval tag
  Tag(intervals_in, num_fixed_ele, fixed_sur_pt_xyz, fixed_sur_ien,
    num_rotated_ele, rotated_sur_pt_xyz, rotated_sur_ien);
}

void Interface_pair::Check_interval(const std::vector<double> &intervals)
{
  SYS_T::print_fatal_if(VEC_T::get_size(intervals) < 2,
    "Error, Interface_pair: illegal 'interval' vector. (1)\n");

  const int n_interval {VEC_T::get_size(intervals) - 1};

  for(int ii{0}; ii<n_interval; ++ii)
  {
    SYS_T::print_fatal_if(intervals[ii] >= intervals[ii + 1],
      "Error, Interface_pair: illegal 'interval' vector. (2)\n");
  }

  if(interface_type == 1)
  {
    SYS_T::print_fatal_if(std::abs(intervals[0] - 0.0) >= 1e-9,
      "Error, Interface_pair: illegal 'interval' vector. (3)");
  }
}

void Interface_pair::Tag(const std::vector<double> &intervals, const int &num_fixed_ele,
  const std::vector<double> &fixed_sur_pt_xyz, const std::vector<int> &fixed_sur_ien,
  const int &num_rotated_ele, const std::vector<double> &rotated_sur_pt_xyz,
  const std::vector<int> &rotated_sur_ien)
{
  fixed_interval_tag.resize(num_fixed_ele);

  for(int ee=0; ee<num_fixed_ele; ++ee)
  {
    Vector_3 ele_centroid(0.0, 0.0, 0.0);
    for(int ii=0; ii<s_nLocBas; ++ii)
    {
      int node = fixed_sur_ien[ee * s_nLocBas + ii];
      ele_centroid.x() += fixed_sur_pt_xyz[3 * node];
      ele_centroid.y() += fixed_sur_pt_xyz[3 * node + 1];
      ele_centroid.z() += fixed_sur_pt_xyz[3 * node + 2];
    }

    ele_centroid *= (1.0 / (double) s_nLocBas);

    fixed_interval_tag[ee] = Group(ele_centroid, intervals);
  }

  rotated_interval_tag.resize(num_rotated_ele);

  for(int ee=0; ee<num_rotated_ele; ++ee)
  {
    Vector_3 ele_centroid(0.0, 0.0, 0.0);
    for(int ii=0; ii<s_nLocBas; ++ii)
    {
      int node = rotated_sur_ien[ee * s_nLocBas + ii];
      ele_centroid.x() += rotated_sur_pt_xyz[3 * node];
      ele_centroid.y() += rotated_sur_pt_xyz[3 * node + 1];
      ele_centroid.z() += rotated_sur_pt_xyz[3 * node + 2];
    }

    ele_centroid *= (1.0 / (double) s_nLocBas);

    rotated_interval_tag[ee] = Group(ele_centroid, intervals);
  }
}

int Interface_pair::Group(const Vector_3 &ele_centroid, const std::vector<double> &intervals)
{
  const int n_interval {VEC_T::get_size(intervals) - 1};

  switch (interface_type)
  {
    case 0:
    {
      for(int ii{0}; ii<n_interval; ++ii)
      {
        if(ele_centroid(T0_axial_direction)>=intervals[ii] && ele_centroid(T0_axial_direction)<=intervals[ii+1])
          return ii;

        SYS_T::print_fatal_if(ii==n_interval && ele_centroid(T0_axial_direction)>intervals[ii+1],
          "Error, Interface_pair: need better intervals end.\n");
      }
    }
    break;

    case 1:
    {
      const double dist {(ele_centroid - T1_surface_centroid).norm2()};
      for(int ii{0}; ii<n_interval; ++ii)
      {
        if(dist>=intervals[ii] && dist<=intervals[ii+1])
          return ii;

        SYS_T::print_fatal_if(ii==n_interval && dist>intervals[ii+1],
          "Error, Interface_pair: need better intervals end.\n");
      }
    }
    break;
    
    default:
      {
      SYS_T::print_fatal("Error, Interface_pair: wrong interface type.\n");
      return -1;
      }
      break;
  }
  return -1;
}

// EOF
