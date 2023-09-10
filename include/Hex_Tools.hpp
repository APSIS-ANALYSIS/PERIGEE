#ifndef HEX_TOOLS_HPP
#define HEX_TOOLS_HPP
#include "DataVecStr.hpp"
#include "Vec_Tools.hpp"
#include "VTK_Tools.hpp"
#include "Math_Tools.hpp"
#include "IIEN.hpp"

#include <vtkQuad.h>
#include <vtkBiQuadraticQuad.h>
#include <vtkHexahedron.h>
#include <vtkTriQuadraticHexahedron.h>

namespace HEX_T
{
  void gen_hex_grid( vtkUnstructuredGrid * const &grid_w,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt,
      const std::vector<int> &ien_array );

  void write_hex_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata, const bool &isXML = true );

void gen_quadrangle_grid( vtkPolyData * const &grid_w,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt,
      const std::vector<int> &ien_array );

void write_quadrangle_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata );

void gen_quadratic_quadrangle_grid( vtkUnstructuredGrid * const &grid_w,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt,
      const std::vector<int> &ien_array );

void write_quadratic_quadrangle_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata );

Vector_3 get_out_normal( const std::string &file,
      const std::vector<double> &vol_ctrlPts,
      const IIEN * const &vol_ien );

  double get_aspect_ratio( const std::vector<double> &coors );

  // ----------------------------------------------------------------
  // ! rest_node: convert the node arrangment of a quadratic hexahedron,
  //              i.e. the local ien, from Gmsh style to VTK style.
  // ----------------------------------------------------------------
  std::vector<int> reset_node (const std::vector<int> &local_ien_msh);
  
}

#endif