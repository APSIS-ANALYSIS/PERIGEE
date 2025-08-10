#ifndef WEDGE_PYRIMID_TOOLS_HPP
#define WEDGE_PYRIMID_TOOLS_HPP

#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Vec_Tools.hpp"

namespace WP_T
{
  Vector_3 get_coor(const int &node_id, const std::vector<double> &ctrlPts);

  bool are_four_points_coplanar(const std::vector<int> &nodes,
    const std::vector<double> &ctrlPts, double &volume);

  // For wedge (prism)  
  std::vector<std::vector<int>> wedge_find_quad_faces(
    const std::vector<int>& node_ids, const std::vector<double>& ctrlPts);

  std::vector<std::vector<int>> wedge_find_side_edges(const std::vector<std::vector<int>> &quad_faces);

  std::vector<std::vector<int>> quad_find_triangule_edges(
    const std::vector<int> &quad_face, const std::vector<int> &side_edge1, const std::vector<int> &side_edge2,
    const std::vector<double> &ctrlPts);

  std::pair<std::vector<int>, std::vector<int>> wedge_find_triangular_faces(
    const std::vector<int> &wedge_nodes, const std::vector<double> &ctrlPts,
    const std::vector<std::vector<int>> &quad_faces, std::vector<std::vector<int>> &side_edges);

  std::vector<std::array<int, 4>> divide_wedge_to_tet(
    const std::vector<int>& node_ids, const std::vector<double>& ctrlPts);

  // For pyrimid
  std::vector<int> pyrimid_find_quad_face(
    const std::vector<int>& node_ids, const std::vector<double>& ctrlPts);

  std::vector<int> quad_find_diagonal(
    const std::vector<int> &quad_face, const std::vector<double> &ctrlPts);

  std::vector<std::array<int, 4>> divide_pyrimid_to_tet(
    const std::vector<int>& node_ids, const std::vector<double>& ctrlPts);
}

#endif