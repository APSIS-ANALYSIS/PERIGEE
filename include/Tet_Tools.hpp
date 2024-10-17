#ifndef TET_TOOLS_HPP
#define TET_TOOLS_HPP
// ==================================================================
// Tet_Tools.hpp
//
// This is a suite of basic tetrahedral element tools with IO and
// basic mesh quality evaluation.
// ==================================================================
#include "DataVecStr.hpp"
#include "Vec_Tools.hpp"
#include "VTK_Tools.hpp"
#include "FE_Tools.hpp"
#include "IIEN.hpp"

#include "vtkTriangle.h"
#include "vtkQuadraticTriangle.h"
#include "vtkTetra.h"
#include "vtkQuadraticTetra.h"

#include "tetgen.h"

namespace TET_T
{
  // ================================================================
  // ===> 1. This set of tools WRITE volumetric mesh of Tetrahedron or
  //         triangle to .vtu file and surface mesh to .vtp file.  
  // ================================================================
  // ----------------------------------------------------------------
  // ! gen_tet_grid: generate the volumetric mesh described by tet
  //                 elements, and pass the data to vtkUnstructuredGrid.
  //   Input: \para grid_w : the vtk object taking the mesh info
  //          \para numpts : the number of grid points
  //          \para numcels : the number of tetrahedral elements
  //          \para pt: xyz coordinates of the linear or quadratic tets,
  //                    length 3 x numpts
  //          \para ien_array : connectivity array, 
  //                            length 4 or 10 x numcels
  // ----------------------------------------------------------------
  void gen_tet_grid( vtkUnstructuredGrid * const &grid_w,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt,
      const std::vector<int> &ien_array );

  // ----------------------------------------------------------------
  // ! write_tet_grid: write the volumetric mesh described by linear 
  //                   or quadratic tetrahedral elements. The routine
  //                   will detect the element type based on the length
  //                   of the ien_array and the numcels.
  //   Input: \para filename : the filename.vtu is the file to be written.
  //          \para numpts : the number of grid points
  //          \para numcels : the number of tetrahedral elements
  //          \para pt: xyz coordinates of the linear or quadratic tets,
  //                    length 3 x numpts
  //          \para ien_array : connectivity array, 
  //                            length 4 or 10 x numcels
  //          \para IOdata : the integer data to be written on cells or nodes 
  //          \para isXML : the flag indicate vtk/vtu format
  // ----------------------------------------------------------------
  void write_tet_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata, const bool &isXML = true );
  
  // ----------------------------------------------------------------
  // ! gen_triangle_grid: generate the surface mesh described by linear triangular
  //                      elements, and pass the data to vtkPolyData.
  //   All input parameters are the same as above, but grid_w is to be
  //   modified by reference. This is convenient for adding additional
  //   fields to the polydata before writing.
  // ----------------------------------------------------------------
  void gen_triangle_grid( vtkPolyData * const &grid_w,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt,
      const std::vector<int> &ien_array );
 
  // ----------------------------------------------------------------
  // ! write_triangle_grid: write the surface mesh described by linear triangular
  //                        elements.
  //   Input: \para filename : the filename.vtp is the file to be written.
  //          \para numpts : the number of grid points
  //          \para numcels : the number of linear triangular elements
  //          \para pt: xyz coordinates of the linear triangles, length 3 numpts
  //          \para ien_array : connectivity array, length 3 numcels
  //          \para IOdata : the integer data to be written on cells or nodes 
  // ----------------------------------------------------------------
  void write_triangle_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata );

  // ----------------------------------------------------------------
  // ! gen_quadratic_triangle_grid: generate the surface mesh described by
  //                                quadratic triangular elements and pass 
  //                                the data to vtkUnstructuredGrid data.
  //   All input parameters are the same as above, but grid_w is to be
  //   modified by reference. This is convenient for adding additional
  //   fields to the unstructured grid before writing.
  // ----------------------------------------------------------------
  void gen_quadratic_triangle_grid( vtkUnstructuredGrid * const &grid_w,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt,
      const std::vector<int> &ien_array );
  
  // ----------------------------------------------------------------
  // ! write_quadratic_triangle_grid: write the surface mesh described 
  //                                  by quadratic triangle elements.
  //   Input: \para filename : the filename.vtu is the file to be written.
  //          \para numpts : the number of grid points
  //          \para numcels : the number of tetrahedral elements
  //          \para pt: xyz coordinates of the quadratic triangles, length 3 numpts
  //          \para ien_array : connectivity array, length 6 numcels
  //          \para IOdata : the integer data to be written on cells or nodes 
  // ----------------------------------------------------------------
  void write_quadratic_triangle_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata );

  // ================================================================
  // 2. Mesh quality measures
  // This set of tools give various measure of tetrahedral element mesh.
  // ================================================================
  // ! get_aspect_ratio:
  //   Input: \para x_i, y_i, z_i for node i=0,1,2,3.
  //          x0 y0 z0 x1 y1 z1 x2 y2 z2 x3 y3 z3
  //   Output: The ratio between the longest edge and the shortest edge.
  //           l_max / l_min.
  //   Lower aspect ratio implies better shape. If the input contains
  //   more than 4 points' coordinate, (e.g. you give the coor of 10-node
  //   tet), the routine will still read the first four nodes. This means
  //   for quadratic tets, the aspect ratio is calculated by extracting
  //   its four vertices and ignoring the mid-edge points.
  // ----------------------------------------------------------------
  double get_aspect_ratio( const std::array<Vector_3, 4> &coors );

  // ----------------------------------------------------------------
  // ! get_out_normal:
  //   This function obtains the unit outward normal vector for a 
  //   triangle vtp or vtu file, with the assumption that the file
  //   describes a flat surface.
  //   The outward normal is calculated based on the first triangle 
  //   element in the vtp file, edge 0-1 and edge 0-2. Making cross
  //   product and compare with the 0-out-of-surface node. If the normal
  //   vector is in the same direction with the 0-out-of-surface line,
  //   correct its direction by multiplying -1.
  //   Note: This function requires one obtain the element index for
  //         the outlet vtp file.
  //   Input: \para file : surface file
  //          \para vol_ctrlPts, the volume mesh control points
  //          \para vol_ien, the volume mesh IEN array
  //   Output: outVec, the outward normal vector
  // ----------------------------------------------------------------
  Vector_3 get_out_normal( const std::string &file,
      const std::vector<double> &vol_ctrlPts,
      const IIEN * const &vol_ien );

  Vector_3 get_out_normal( const std::string &file,
      const std::vector<double> &vol_ctrlPts,
      const IIEN * const &vol_ien,
      const double &num_shifted_node,
      const double &num_shifted_elem );
  
  // ================================================================
  // 3. TetGen interface
  // ----------------------------------------------------------------
  // This set of tools convert the tetgenio object to vtu grid and vtp
  // surface files with markers for the preprocess code to read.
  // ================================================================
  // Input: \para meshout: tetgenio object containing the mesh points 
  //                      and cell connectivity.
  //        \para fName: fName.vtu will be write on disk as the volumetric
  //                     grid.
  // ----------------------------------------------------------------
  void tetgenio2vtu( const tetgenio &meshout, const std::string &fName );

  // ----------------------------------------------------------------  
  // Input: \para meshout: tetgenio object containing the mesh points 
  //                      and cell connectivity.
  //        \para fName: fName.vtp will be write on disk as the surface 
  //                     grid.
  //        \para bcmarker: the marker input in the tetgen input file
  //                        specifying the boundary domain.
  // Output: fName.bcmarker.vtp will be write on disk for the boundary
  //         with marker bcmarker.
  // ----------------------------------------------------------------
  void tetgenio2vtp( const tetgenio &meshout, const std::string &fName,
      const int &bcmarker );

  // ================================================================
  // 4. Tet4 class defines a four node linear tetrahedron object with 
  // basic manipulations including checking the size, the aspect ratio,
  // and check the face index for given three nodal indices.
  // The object is defined by 12 double data: 
  //             x0, y0, z0, x1, ..., z3
  // ================================================================
  class Tet4
  {
    public:
      // Default constructor: generate a default tetrahedron in the
      // following form
      //          x y z
      // node 0 : 0 0 0
      // node 1 : 1 0 0
      // node 2 : 0 1 0
      // node 3 : 0 0 1
      Tet4();

      // Generate a tetrahedron with the x-y-z coordinates of the
      // four points given by nodes-vector
      Tet4( const std::array<Vector_3, 4> &in_nodes );

      Tet4( const std::array<Vector_3, 4> &input_pts, 
          const int &ien0, const int &ien1, 
          const int &ien2, const int &ien3 );

      virtual ~Tet4();

      // Changes the gindex array only. It is used to determine
      // the face index, e.g. in ElemBC_3D::resetSurIEN_outwardnormal
      // function.
      void reset( const int &ien0, const int &ien1,
          const int &ien2, const int &ien3 );

      void reset( const std::array<int,4> &ienval )
      { reset(ienval[0], ienval[1], ienval[2], ienval[3]); }

      // Reset the object content includign the global index and the point
      // coordindates
      void reset( const std::vector<double> &ctrlPts,
          const IIEN * const &ien_pt, const int &ee );

      double get_aspect_ratio() const;

      double get_volume() const;

      // Get the circumscribing sphere's DIAMETER.
      // This is a measure of the mesh size
      double get_diameter() const;

      // Given the face node IEN indices, determine the face id is determined by
      // its opposite node number.
      int get_face_id(const int &n0, const int &n1, const int &n2) const;

      int get_face_id( const std::array<int,3> &nn ) const
      { return get_face_id(nn[0], nn[1], nn[2]); }

      void print_info() const;

    private:
      std::array<Vector_3, 4> pts;

      int gindex[4];
  };

  // ================================================================
  // 5. Tetrahedral mesh checker. 
  //    This routine will read in the mesh information: control points
  //    and IEN array, and check each element's quality and print
  //    necessary information.
  //    cpts: list of control points;
  //    ienptr: IEN array
  //    nelem: total number of element
  //    crit_aspect_ratio: the element above this value will be counted.
  // ================================================================
  void tetmesh_check(const std::vector<double> &cpts,
      const IIEN * const &ienptr, const int &nelem,
      const double &crit_aspect_ratio = 3.5 );
}

#endif
