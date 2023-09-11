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

  class Hex8
  {
    public:
      // Default constructor: generate a default hexahedron in the
      // following form
      //           x  y  z
      // node 0 : -1 -1 -1
      // node 1 :  1 -1 -1
      // node 2 :  1  1 -1
      // node 3 : -1  1 -1
      // node 4 : -1 -1  1
      // node 5 :  1 -1  1
      // node 6 :  1  1  1
      // node 7 : -1  1  1
      Hex8();

      // Generate a hexahedron with the x-y-z coordinates of the
      // eight points given by nodes-vector
      Hex8( const std::vector<double> &in_nodes );

      virtual ~Hex8();

      void reset( const std::vector<double> &nodes );

      private:
      double pts[24];

      int gindex[8];
  };
}

#endif