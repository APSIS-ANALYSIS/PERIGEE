#ifndef HEX_TOOLS_HPP
#define HEX_TOOLS_HPP
// ============================================================================
// Hex_Tools.hpp
//
// This is a suite of basic hexahedral element tools with IO and basic mesh 
// quality evaluation.

// We use 27-node hexahedron and 9-node quadrilateral, but NOT 20-node hexahedron
// or 8-node quadrilateral, as the quadratic
// elements at present.
// ============================================================================
#include "Tet_Tools.hpp"

#include <vtkQuad.h>
#include <vtkBiQuadraticQuad.h>
#include <vtkHexahedron.h>
#include <vtkTriQuadraticHexahedron.h>

namespace HEX_T
{
  // ==========================================================================
  // ===> 1. This set of tools WRITE volumetric mesh of Hexahedron or
  //         quadrilateral to .vtu file and surface mesh to .vtp file.  
  // ==========================================================================
  // --------------------------------------------------------------------------
  // ! gen_hex_grid: generate the volumetric mesh described by hex
  //                 elements, and pass the data to vtkUnstructuredGrid.
  //   Input: \para grid_w : the vtk object taking the mesh info
  //          \para numpts : the number of grid points
  //          \para numcels : the number of hexahedral elements
  //          \para pt: xyz coordinates of the linear or quadratic hexes,
  //                    length 3 x numpts
  //          \para ien_array : connectivity array, length is 8 or 27 x numcels
  // --------------------------------------------------------------------------
  void gen_hex_grid( vtkUnstructuredGrid * const &grid_w,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt,
      const std::vector<int> &ien_array );

  // --------------------------------------------------------------------------
  // ! write_tet_grid: write the volumetric mesh described by linear 
  //                   or quadratic hexahedral elements. The routine
  //                   will detect the element type based on the length
  //                   of the ien_array and the numcels.
  //   Input: \para filename : the filename.vtu is the file to be written.
  //          \para numpts : the number of grid points
  //          \para numcels : the number of hexahedral elements
  //          \para pt: xyz coordinates of the linear tets, length 3 x numpts
  //          \para ien_array : connectivity array, length 8 or 27 x numcels
  //          \para IOdata : the integer data to be written on cells or nodes 
  //          \para isXML : the flag indicate vtk/vtu format
  // --------------------------------------------------------------------------
  void write_hex_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata, const bool &isXML = true );

  // --------------------------------------------------------------------------
  // ! gen_quad_grid: generate the surface mesh described by linear quadrilateral 
  //                  elements, and pass the data to vtkPolyData.
  //   All input parameters are the same as above, but grid_w is to be
  //   modified by reference. This is convenient for adding additional
  //   fields to the polydata before writing.
  // --------------------------------------------------------------------------
  void gen_quad_grid( vtkPolyData * const &grid_w,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt,
      const std::vector<int> &ien_array );

  // --------------------------------------------------------------------------
  // ! write_quad_grid: write the surface mesh described by linear quadrilateral
  //                    elements.
  //   Input: \para filename : the filename.vtp is the file to be written.
  //          \para numpts : the number of grid points
  //          \para numcels : the number of linear quadrilateral elements
  //          \para pt: xyz coordinates of the linear quads, length 3 numpts
  //          \para ien_array : connectivity array, length 4 numcels
  //          \para IOdata : the integer data to be written on cells or nodes 
  // --------------------------------------------------------------------------
  void write_quad_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata );

  // --------------------------------------------------------------------------
  // ! gen_quadratic_quad_grid: generate the surface mesh described by
  //                            quadratic quadrilateral elements and pass 
  //                            the data to vtkUnstructuredGrid data.
  //   All input parameters are the same as above, but grid_w is to be
  //   modified by reference. This is convenient for adding additional
  //   fields to the unstructured grid before writing.
  // --------------------------------------------------------------------------
  void gen_quadratic_quad_grid( vtkUnstructuredGrid * const &grid_w,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt,
      const std::vector<int> &ien_array );

  // --------------------------------------------------------------------------
  // ! write_quadratic_quad_grid: write the surface mesh described by
  //                              quadratic quadrilateral elements.
  //   Input: \para filename : the filename.vtu is the file to be written.
  //          \para numpts : the number of grid points
  //          \para numcels : the number of quadrilateral elements
  //          \para pt: xyz coordinates of the quadratic quads, length 3 numpts
  //          \para ien_array : connectivity array, length 9 numcels
  //          \para IOdata : the integer data to be written on cells or nodes 
  // --------------------------------------------------------------------------
  void write_quadratic_quad_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata );

  // ==========================================================================
  // 2. Mesh quality measures
  // This set of tools give various measure of hexahedral element mesh.
  // ==========================================================================
  // ! get_out_normal:
  //   This function obtains the unit outward normal vector for a 
  //   quadrilateral vtp or vtu file, with the assumption that the file
  //   describes a flat surface.
  //   The outward normal is calculated based on the first quadrilateral 
  //   element in the vtp file, edge 0-1 and edge 0-3. Making cross
  //   product and compare with the inward vector, which starts at the
  //   center of the quadrilater and points to the center of hexaheron where
  //   the quad is attached. If the inner product of the normal vector
  //   and the inward vector, correct the direction of normal vector by 
  //   multiplying -1.
  //   Note: This function requires one obtain the element index for
  //         the outlet vtp file.
  //   Input: \para file : surface file
  //          \para vol_ctrlPts, the volume mesh control points
  //          \para vol_ien, the volume mesh IEN array
  //   Output: outVec, the outward normal vector
  // --------------------------------------------------------------------------
  Vector_3 get_out_normal( const std::string &file,
      const std::vector<double> &vol_ctrlPts,
      const IIEN * const &vol_ien );

  // --------------------------------------------------------------------------
  // ! get_aspect_ratio:
  //   Input: \para x_i, y_i, z_i for node i=0,1,2,3,4,5,6,7.
  //          x0 y0 z0 x1 y1 z1 x2 y2 z2 x3 y3 z3 ... x7 y7 z7
  //   Output: The ratio between the longest edge and the shortest edge.
  //           l_max / l_min.
  //   Lower aspect ratio implies better shape. If the input contains
  //   more than 8 points' coordinate, (e.g. you give the coor of 27-node
  //   hex), the routine will still read the first eight nodes. This means
  //   for quadratic hexes, the aspect ratio is calculated by extracting
  //   its eight vertices and ignoring the mid-edge points.
  // --------------------------------------------------------------------------
  double get_aspect_ratio( const std::array<Vector_3, 8> &pt );

  // ==========================================================================
  // 3. Hex8 class defines a four node linear tetrahedron object with 
  // basic manipulations including checking the size, the aspect ratio.
  // The object is defined by 12 double data: 
  //             x0, y0, z0, x1, ..., z7
  // ==========================================================================
  class Hex8
  {
    public:
      // ----------------------------------------------------------------------
      // Default constructor: generate a default hexahedron in the following form
      //           x  y  z
      // node 0 : -1 -1 -1
      // node 1 :  1 -1 -1
      // node 2 :  1  1 -1
      // node 3 : -1  1 -1
      // node 4 : -1 -1  1
      // node 5 :  1 -1  1
      // node 6 :  1  1  1
      // node 7 : -1  1  1
      // ----------------------------------------------------------------------
      Hex8();

      // Generate a hexahedron with the x-y-z coordinates of the
      // eight points given by nodes-vector
      Hex8( const std::array<Vector_3, 8> &in_nodes );

      Hex8( const std::array<Vector_3, 8> &input_pts,
          const int &ien0, const int &ien1, const int &ien2,
          const int &ien3, const int &ien4, const int &ien5,
          const int &ien6, const int &ien7 );

      virtual ~Hex8();
      
      // Changes the gindex array only. It is used to determine the face index, 
      // e.g. in ElemBC_3D::resetSurIEN_outwardnormal function.
      void reset( const int &ien0, const int &ien1,
          const int &ien2, const int &ien3, const int &ien4,
          const int &ien5, const int &ien6, const int &ien7 );

      void reset( const std::vector<double> &ctrlPts,
          const IIEN * const &ien_ptr, const int &ee );

      void print_info() const;

      double get_aspect_ratio() const;

      // This is a measure of the mesh size, defined by
      // the max length of four body diagonals
      double get_diameter() const;

      double get_volume() const;

      // Given the face node IEN indices, determine the face id
      int get_face_id(const int &n0, const int &n1, const int &n2, const int &n3) const;

    private:
      std::array<Vector_3, 8> pts;

      int gindex[8];
  };

  // ================================================================
  // 4. Hexahedral mesh checker. 
  //    This routine will read in the mesh information: control points
  //    and IEN array, and check each element's quality and print
  //    necessary information.
  //    cpts: list of control points;
  //    ienptr: IEN array
  //    nelem: total number of element
  //    crit_aspect_ratio: the element above this value will be counted.
  // ================================================================
  void hexmesh_check(const std::vector<double> &cpts,
      const IIEN * const &ienptr, const int &nelem,
      const double &crit_aspect_ratio = 3.5 );
}

#endif
