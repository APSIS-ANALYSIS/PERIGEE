#include "Wedge_Pyrimid_Tools.hpp"

Vector_3 WP_T::get_coor(const int &node_id, const std::vector<double> &ctrlPts)
{
  return Vector_3( ctrlPts[3*node_id], ctrlPts[3*node_id+1], ctrlPts[3*node_id+2] );
}

bool WP_T::are_four_points_coplanar(const std::vector<int> &nodes,
  const std::vector<double> &ctrlPts, double &volume)
{
  Vector_3 p1 = WP_T::get_coor(nodes[0], ctrlPts);
  Vector_3 p2 = WP_T::get_coor(nodes[1], ctrlPts);
  Vector_3 p3 = WP_T::get_coor(nodes[2], ctrlPts);
  Vector_3 p4 = WP_T::get_coor(nodes[3], ctrlPts);

  Vector_3 v12 = p2 - p1;
  Vector_3 v13 = p3 - p1;
  Vector_3 v14 = p4 - p1;

  Vector_3 cross_12_13 = Vec3::cross_product(v12, v13);

  volume = std::abs( cross_12_13.dot_product(v14) );

  double min_mag = VEC_T::min(std::vector<double> {v12.norm2(), v13.norm2(), v14.norm2()});

  return volume < 0.333 * min_mag * min_mag * min_mag;
}

std::vector<std::vector<int>> WP_T::wedge_find_quad_faces(
const std::vector<int>& node_ids, const std::vector<double>& ctrlPts)
{
  std::vector<std::vector<int>> quad_faces (3, std::vector<int>(4, -1));

  std::vector<std::vector<int>> all_cases {};
  std::vector<double> volumes {};
  double vol = 0.0;

  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=ii+1; jj<4; ++jj)
    {
      for(int kk=jj+1; kk<5; ++kk)
      {
        for(int ll=kk+1; ll<6; ++ll)
        {
          std::vector<int> four_nodes = {node_ids[ii], node_ids[jj], node_ids[kk], node_ids[ll]};
          WP_T::are_four_points_coplanar(four_nodes, ctrlPts, vol);

          all_cases.push_back(four_nodes);
          volumes.push_back(vol);
        }
      }
    }
  }

  double min_volume = VEC_T::min(volumes);
  int pos = VEC_T::get_pos(volumes, min_volume);

  quad_faces[0] = all_cases[pos];
  VEC_T::erase_pos(all_cases, pos);
  VEC_T::erase_pos(volumes, pos);

  min_volume = VEC_T::min(volumes);
  pos = VEC_T::get_pos(volumes, min_volume);

  quad_faces[1] = all_cases[pos];
  VEC_T::erase_pos(all_cases, pos);
  VEC_T::erase_pos(volumes, pos);

  min_volume = VEC_T::min(volumes);
  pos = VEC_T::get_pos(volumes, min_volume);

  quad_faces[2] = all_cases[pos];
  VEC_T::erase_pos(all_cases, pos);
  VEC_T::erase_pos(volumes, pos);

  return quad_faces;
}

std::vector<std::vector<int>> WP_T::wedge_find_side_edges(const std::vector<std::vector<int>> &quad_faces)
{
  std::vector<std::vector<int>> side_edges;

  for(int ii=0; ii<2; ++ii)
  {
    for(int jj=ii+1; jj<3; ++jj)
    {
      std::vector<int> face1 = quad_faces[ii];
      std::vector<int> face2 = quad_faces[jj];

      std::vector<int> common_edge = VEC_T::intersection(face1, face2);

      side_edges.push_back(common_edge);

      if(side_edges.size() == 3) // A wedge has exactly 3 side edges
        return side_edges;
    }
  }

  if(side_edges.size() != 3)
    throw std::runtime_error("Error: Expected exactly 3 side edges for a wedge.");
        
  return side_edges;
}

std::vector<std::vector<int>> WP_T::quad_find_triangule_edges(
const std::vector<int> &quad_face, const std::vector<int> &side_edge1, const std::vector<int> &side_edge2,
const std::vector<double> &ctrlPts)
{
  std::vector<int> quad = quad_face;
  VEC_T::sort_unique_resize(quad);

  std::vector<int> temp_quad = side_edge1;
  VEC_T::insert_end(temp_quad, side_edge2);
  VEC_T::sort_unique_resize(temp_quad);

  std::vector<int> temp_empty = VEC_T::set_diff(quad, side_edge1);
  temp_empty = VEC_T::set_diff(temp_empty, side_edge2);

  if(VEC_T::is_equal(quad, temp_quad) && temp_empty.size() == 0)
  {
    Vector_3 v1 = WP_T::get_coor(side_edge1[1], ctrlPts) - WP_T::get_coor(side_edge1[0], ctrlPts);
    Vector_3 v2 = WP_T::get_coor(side_edge2[1], ctrlPts) - WP_T::get_coor(side_edge2[0], ctrlPts);

    if(v1.dot_product(v2) > 0.0)
    {
      std::vector<int> tri_edge1 = {side_edge1[0], side_edge2[0]};
      std::vector<int> tri_edge2 = {side_edge1[1], side_edge2[1]};

      return std::vector<std::vector<int>>{tri_edge1, tri_edge2};
    }
    else
    {
      std::vector<int> tri_edge1 = {side_edge1[0], side_edge2[1]};
      std::vector<int> tri_edge2 = {side_edge1[1], side_edge2[0]};

      return std::vector<std::vector<int>>{tri_edge1, tri_edge2};
    }  
  }
  else
    return {};
}

std::pair<std::vector<int>, std::vector<int>> WP_T::wedge_find_triangular_faces(
const std::vector<int> &wedge_nodes, const std::vector<double> &ctrlPts,
const std::vector<std::vector<int>> &quad_faces, std::vector<std::vector<int>> &side_edges)
{
  // std::vector<int> tri_face1 {};
  std::vector<std::vector<int>> triangle_edges {};

  for(const auto &quad_face : quad_faces)
  {
    for(int ii=0; ii<2; ++ii)
    {
      for(int jj=ii+1; jj<3; ++jj)
      {
        std::vector<int> side_edge1 = side_edges[ii];
        std::vector<int> side_edge2 = side_edges[jj];
        
        auto tri_edge = WP_T::quad_find_triangule_edges(quad_face, side_edge1, side_edge2, ctrlPts);

        if(VEC_T::get_size(tri_edge) != 0)
        {
          for(const auto &edge : tri_edge)
            triangle_edges.push_back(edge);
        }
      }
    }
  }

  // std::cout<< "Triangle edges found: " << triangle_edges.size() << std::endl;
  // for(const auto &edge : triangle_edges)
  // {
  //   VEC_T::print(edge);
  // }

  if(triangle_edges.size() != 6)
    throw std::runtime_error("Error: Expected exactly 6 triangular edges for a wedge.");

  // for(int ii=1; ii<6; ++ii)
  // {
  //   std::vector<int> inter_nodes = VEC_T::intersection(triangle_edges[0], triangle_edges[ii]);
  //   if(VEC_T::get_size(inter_nodes) != 0)
  //   {
  //     tri_face1 = triangle_edges[0]; 
  //     VEC_T::insert_end(tri_face1, triangle_edges[ii]);
  //     VEC_T::sort_unique_resize(tri_face1);
  //     break;
  //   }
  // }

  // std::vector<int> tri_face2 = VEC_T::set_diff(wedge_nodes, tri_face1);

  Vector_3 ve1 = WP_T::get_coor(side_edges[0][1], ctrlPts) - WP_T::get_coor(side_edges[0][0], ctrlPts);
  Vector_3 ve2 = WP_T::get_coor(side_edges[1][1], ctrlPts) - WP_T::get_coor(side_edges[1][0], ctrlPts);
  Vector_3 ve3 = WP_T::get_coor(side_edges[2][1], ctrlPts) - WP_T::get_coor(side_edges[2][0], ctrlPts);

  int temp = -1;
  if(ve1.dot_product(ve2) < 0.0)
  {
    temp = side_edges[1][0];
    side_edges[1][0] = side_edges[1][1];
    side_edges[1][1] = temp;
  }

  if(ve1.dot_product(ve3) < 0.0)
  {
    temp = side_edges[2][0];
    side_edges[2][0] = side_edges[2][1];
    side_edges[2][1] = temp;
  }
  
  std::vector<int> tri_face1 = {side_edges[0][0], side_edges[1][0], side_edges[2][0]};
  std::vector<int> tri_face2 = {side_edges[0][1], side_edges[1][1], side_edges[2][1]};

  return std::make_pair(tri_face1, tri_face2);
}

std::vector<std::array<int, 4>> WP_T::divide_wedge_to_tet(
const std::vector<int>& node_ids, const std::vector<double>& ctrlPts)
{
  SYS_T::print_fatal_if(VEC_T::get_size(node_ids) != 6,
    "Error: Expected exactly 6 node IDs for a wedge.");

  auto quad_faces = WP_T::wedge_find_quad_faces(node_ids, ctrlPts);

  // std::cout << "Quadrilateral faces found: " << quad_faces.size() << std::endl;
  // for(const auto &face : quad_faces)
  // {
  //   VEC_T::print(face);
  // }

  std::vector<std::vector<int>> side_edges;
  try
  {
    side_edges = WP_T::wedge_find_side_edges(quad_faces);
  }
  catch(const std::exception& e)
  {
    // std::cerr << e.what() << '\n';
    throw ;
  }

  // std::cout << "Side edges found: " << side_edges.size() << std::endl;
  // for(const auto &side_edge : side_edges)
  // {
  //   VEC_T::print(side_edge);
  // }

  std::pair<std::vector<int>, std::vector<int>> triangular_faces;
  try
  {
    triangular_faces = WP_T::wedge_find_triangular_faces(node_ids, ctrlPts, quad_faces, side_edges);
  }
  catch(const std::exception& e)
  {
    // std::cerr << e.what() << '\n';
    throw ;
  }
  
  std::vector<int> face1 = triangular_faces.first;
  std::vector<int> face2 = triangular_faces.second;

  std::vector<std::array<int, 4>> tets;

  tets.push_back({face1[0], face1[1], face1[2], face2[0]});

  tets.push_back({face1[1], face1[2], face2[0], face2[1]});

  tets.push_back({face1[2], face2[0], face2[1], face2[2]});

  return tets;
}

std::vector<int> WP_T::pyrimid_find_quad_face(
const std::vector<int>& node_ids, const std::vector<double>& ctrlPts)
{
  std::vector<std::vector<int>> all_cases {};
  std::vector<double> volumes {};
  double vol = 0.0;

  for(int ii=0; ii<2; ++ii)
  {
    for(int jj=ii+1; jj<3; ++jj)
    {
      for(int kk=jj+1; kk<4; ++kk)
      {
        for(int ll=kk+1; ll<5; ++ll)
        {
          std::vector<int> four_nodes = {node_ids[ii], node_ids[jj], node_ids[kk], node_ids[ll]};
          WP_T::are_four_points_coplanar(four_nodes, ctrlPts, vol);

          all_cases.push_back(four_nodes);
          volumes.push_back(vol);
        }
      }
    }
  }

  double min_volume = VEC_T::min(volumes);
  int pos = VEC_T::get_pos(volumes, min_volume);

  return all_cases[pos];
}

std::vector<int> WP_T::quad_find_diagonal(
const std::vector<int> &quad_face, const std::vector<double> &ctrlPts)
{
  std::vector<int> diagonal {};

  double max_area = 0.0;

  for(int ii=1; ii<4; ++ii)
  {
    std::vector<int> temp_diagonal = {quad_face[0], quad_face[ii]};
    std::vector<int> other_nodes = VEC_T::set_diff(quad_face, temp_diagonal);

    Vector_3 vec_diag = WP_T::get_coor(temp_diagonal[1], ctrlPts) - WP_T::get_coor(temp_diagonal[0], ctrlPts);

    Vector_3 vec1 = WP_T::get_coor(other_nodes[0], ctrlPts) - WP_T::get_coor(temp_diagonal[0], ctrlPts);
    Vector_3 vec2 = WP_T::get_coor(other_nodes[1], ctrlPts) - WP_T::get_coor(temp_diagonal[0], ctrlPts);

    double area = Vec3::cross_product(vec_diag, vec1).norm2() + Vec3::cross_product(vec_diag, vec2).norm2();

    if(area > max_area)
    {
      max_area = area;
      diagonal = temp_diagonal;
    }
  }

  return diagonal;
}

std::vector<std::array<int, 4>> WP_T::divide_pyrimid_to_tet(
const std::vector<int>& node_ids, const std::vector<double>& ctrlPts)
{
  SYS_T::print_fatal_if(VEC_T::get_size(node_ids) != 5,
    "Error: Expected exactly 5 node IDs for a pyrimid.");

  auto quad_face = WP_T::pyrimid_find_quad_face(node_ids, ctrlPts);
  auto diagonal = WP_T::quad_find_diagonal(quad_face, ctrlPts);

  auto peak_node = VEC_T::set_diff(node_ids, quad_face);

  // Peak node and the diagonal nodes
  VEC_T::insert_end(peak_node, diagonal);
  auto other_nodes = VEC_T::set_diff(node_ids, peak_node);

  std::array<int, 4> tet1 = {peak_node[0], diagonal[0], diagonal[1], other_nodes[0]};
  std::array<int, 4> tet2 = {peak_node[0], diagonal[1], diagonal[0], other_nodes[1]};
  
  return std::vector<std::array<int, 4>> {tet1, tet2};
}