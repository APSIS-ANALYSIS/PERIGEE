#ifndef TET_TOOLS_HPP
#define TET_TOOLS_HPP
// ==================================================================
// Tet_Tools.hpp
//
// This is a suite of basic tetrahedral element tools with IO and
// basic mesh quality evaluation.
// 
// Author: Ju Liu, liujuy@gmail.com
// ==================================================================
#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Vector_3.hpp"
#include "IIEN.hpp"

#include "vtkCellLocator.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkTriangle.h"
#include "vtkQuadraticTriangle.h"
#include "vtkTetra.h"
#include "vtkQuadraticTetra.h"
#include "vtkGenericCell.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLGenericDataObjectReader.h"

#include "tetgen.h"

namespace TET_T
{
  // ================================================================
  // ===> 1. The first set of tools READ volumetric mesh from .vtu 
  //         file and surface mesh from .vtp file.  
  // ================================================================
  // --------------------------------------------------------------
  // ! read_vtu_grid: read the mesh generated from other software
  //                  in the .vtu file. The mesh file is assumed to be a
  //                  tetrahedral mesh with 4 or 10 nodes;
  //                  or a triangle mesh with 6 nodes; 
  //                  otherwise, an error message will be thrown.
  //   Input:  \para filename : the filename ending with .vtu
  //   Output: \para numpts: the number of grid points
  //           \para numcels: the number of cells
  //           \para pt: xyz coordinate of the grids, 
  //                     length is 3 x numpts
  //           \para ien_array: the connectivity array, 
  //                            length is 4, 10, or 6 x numcels
  //   Note: reload the function with output in int type instead
  //         of the VTK builtin type vtkIdType by doing a static_cast.
  // ----------------------------------------------------------------
  void read_vtu_grid( const std::string &filename,
      int &numpts, int &numcels,
      std::vector<double> &pt, std::vector<int> &ien_array );

  // ----------------------------------------------------------------
  // read integer / double cell or point data from file in vtu / vtp format.
  // Input: \para filename the vtk file name
  //        \para dataname the data property name
  // Output: \para data the data associated with cell/point
  // ----------------------------------------------------------------
  std::vector<int> read_int_CellData( const std::string &filename, 
      const std::string &dataname );

  std::vector<double> read_double_CellData( const std::string &filename, 
      const std::string &dataname );

  std::vector<int> read_int_PointData( const std::string &filename, 
      const std::string &dataname );

  std::vector<double> read_double_PointData( const std::string &filename, 
      const std::string &dataname );

  // ----------------------------------------------------------------
  // ! read_vtu_grid: read the mesh info just exactly the same way
  //                  as the original read_vtu_grid. In addition,
  //                  this function reads the <physical tag> as an 
  //                  additional output data.
  //   Input:  \para filename : the filename ending with .vtu
  //   Output: \para numpts: the number of grid points
  //           \para numcels: the number of cells
  //           \para pt: xyz coordinate of the grids, length is 3 numpts
  //           \para ien_array: the connectivity array, 
  //                            length is 4, 10, or 6 x numcels
  //           \para phy_tag: the tag of the elements, length is numcels
  //   This function is specifically designed for FSI or multi-domain 
  //   problems, where we need a tag to identify the physical domain.
  // ----------------------------------------------------------------
  void read_vtu_grid( const std::string &filename,
      int &numpts, int &numcels,
      std::vector<double> &pt, std::vector<int> &ien_array,
      std::vector<int> &phy_tag );

  // ----------------------------------------------------------------
  // ! read_vtu_grid: read the surface mesh generated from other software
  //                  in .vtu files. The mesh file is assumed to be a VTK
  //                  triangle grid (type 22 in VTK cell type); 
  //                  otherwise, an error message will be thrown.
  //   Input:  \para filename : the file name ending with .vtu
  //   Output: \para numpts: the number of grid points
  //           \para numcels: the number of triangle cells
  //           \para pt: xyz coordinate of the triangles, 
  //                     length is 3 x numpts.
  //           \para ien_array: the connectivity array, 
  //                            length is 6 x numcels.
  //           \para global_node_index: the mapping from local nodal
  //                                    index to global nodal index
  //           \para global_ele_index: the mapping from the triangle
  //                                   cell to its corresponding tet
  //                                   cell index
  // ----------------------------------------------------------------
  void read_vtu_grid( const std::string &filename,
      int &numpts, int &numcels,
      std::vector<double> &pt, std::vector<int> &ien_array,
      std::vector<int> &global_node_index,
      std::vector<int> &global_elem_index );


  // ----------------------------------------------------------------
  // ! read_vtp_grid: read the surface mesh from a .vtp file. The mesh
  //                  file is assumed to be a VTK trangle grid with 3
  //                  nodes. 
  //                  We do not read in the associated volumetric nodal
  //                  or element indices in this routine. Hence, this 
  //                  is used primarily for reading Dirichlet type BC 
  //                  information.
  //                  
  //   Input:  \para filename : the file name ending with .vtp
  //   Output: \para numpts : the number of grid points
  //           \para numcels : the number of triangle cells
  //           \para pt : the xyz-coordinate of grid points
  //           \para ien_array: the connectivity array for the triangles,
  //                            length is 3 x numcels
  // ----------------------------------------------------------------
  void read_vtp_grid( const std::string &filename,
      int &numpts, int &numcels,
      std::vector<double> &pt, std::vector<int> &ien_array );


  // ----------------------------------------------------------------
  // ! read_vtp_grid: read the surface mesh generated from other software
  //                  in .vtp files. The mesh file is assumed to be a VTK
  //                  triangle grid (type 5 in VTK cell type); 
  //                  otherwise, an error message will be thrown.
  //   Input:  \para filename : the file name ending with .vtp
  //   Output: \para numpts: the number of grid points
  //           \para numcels: the number of triangle cells
  //           \para pt: xyz coordinate of the triangles, 
  //                     length is 3 x numpts.
  //           \para ien_array: the connectivity array, 
  //                            length is 3 x numcels.
  //           \para global_node_index: the mapping from local nodal
  //                                    index to global nodal index
  //           \para global_ele_index: the mapping from the triangle
  //                                   cell to its corresponding tet
  //                                   cell index
  // ----------------------------------------------------------------
  void read_vtp_grid( const std::string &filename,
      int &numpts, int &numcels,
      std::vector<double> &pt, std::vector<int> &ien_array,
      std::vector<int> &global_node_index,
      std::vector<int> &global_ele_index );


  // ================================================================
  // ===> 2. The second set of tools WRITE volumetric mesh to .vtu 
  //         file and surface mesh to .vtp file.  
  // ================================================================
  // ----------------------------------------------------------------
  // ! gen_tet_grid: generate the volumetric mesh described by tet
  //                 elements, and pass the data to vtkUnstructuredGrid.
  //   Input: \para grid_w : the vtk object taking the mesh info
  //          \para numpts : the number of grid points
  //          \para numcels : the number of tetrahedral elements
  //          \para pt: xyz coordinates of the linear tets, length 
  //                    3 x numpts
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
  //                   will detect the element type from the length of
  //                   the ien_array and the numcels.
  //   Input: \para filename : the filename.vtu is the file to be written.
  //          \para numpts : the number of grid points
  //          \para numcels : the number of tetrahedral elements
  //          \para pt: xyz coordinates of the linear tets, length 
  //                    3 x numpts
  //          \para ien_array : connectivity array, 
  //                            length 4 or 10 x numcels
  //          \para node_idx : the nodal index associated with nodes
  //          \para elem_idx : the element index asscoiated with cells
  //
  //   If node_idx, elem_idx are not given, a natural numbering will
  //   be generated for nodes and cells.
  //   Additional argument, 
  //          \para phy_tag : the tag indicating different physical
  //                          domains.
  //          \para isXML : the flag indicate vtk/vtu format
  // ----------------------------------------------------------------
  void write_tet_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, const std::vector<int> &ien_array );

  void write_tet_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, const std::vector<int> &ien_array,
      const std::vector<int> &node_idx, const std::vector<int> &elem_idx );

  void write_tet_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, const std::vector<int> &ien_array,
      const std::vector<int> &node_idx, const std::vector<int> &elem_idx,
      const std::vector<int> &phy_tag, const bool &isXML );

  // ----------------------------------------------------------------
  // ! write_tet_grid : write a volumetric mesh with 
  //                    1) NodalIndex; 2) ElemIndex; 3) Phy_tag
  //   Additional Input compared with original write_tet_grid: 
  //          \para phy_tag : the vector of physical domain tags
  //          \para isXML : the flag that determines write as vtk or
  //                        vtu file
  //          \para start_cell_index : the starting cell/element index
  //                        in this subdomain, by default 0
  // ----------------------------------------------------------------
  void write_tet_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, const std::vector<int> &ien_array,
      const std::vector<int> &phy_tag, const bool &isXML,
      const int &start_cell_index = 0 );


  // ----------------------------------------------------------------
  // ! write_triangle_grid: write the surface mesh described by triangle
  //                        elements.
  //   Input: \para filename : the filename.vtp is the file to be written.
  //          \para numpts : the number of grid points
  //          \para numcels : the number of tetrahedral elements
  //          \para pt: xyz coordinates of the linear tets, length 3 numpts
  //          \para ien_array : connectivity array, length 3 numcels
  //          \para nodal_index : the point data to be written
  //          \para elem_index : the element index to be written
  // ----------------------------------------------------------------
  void write_triangle_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<int> &global_node_index,
      const std::vector<int> &global_ele_index );


  // ----------------------------------------------------------------
  // ! gen_triangle_grid: generate the surface mesh described by triangle
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
  // ! write_quadratic_triangle_grid: write the surface mesh described 
  //                                  by quadratic triangle elements.
  //   Input: \para filename : the filename.vtu is the file to be written.
  //          \para numpts : the number of grid points
  //          \para numcels : the number of tetrahedral elements
  //          \para pt: xyz coordinates of the linear tets, length 3 numpts
  //          \para ien_array : connectivity array, length 3 numcels
  //          \para nodal_index : the point data to be written
  //          \para elem_index : the element index to be written
  // ----------------------------------------------------------------
  void write_quadratic_triangle_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<int> &global_node_index,
      const std::vector<int> &global_ele_index );


  // ----------------------------------------------------------------
  // ! gen_quadratic_triangle_grid: generate the surface mesh described by
  //                                triangle elements and pass the data
  //                                to vtkUnstructuredGrid data.
  //   All input parameters are the same as above, but grid_w is to be
  //   modified by reference. This is convenient for adding additional
  //   fields to the unstructured grid before writing.
  // ----------------------------------------------------------------
  void gen_quadratic_triangle_grid( vtkUnstructuredGrid * const &grid_w,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt,
      const std::vector<int> &ien_array );


  // ----------------------------------------------------------------
  // ! write_triangle_grid: write the surface mesh described by triangle
  //                        elements with two element index arrays.
  //   The input parameter is identical to the former one, with one
  //   more array, global_ele_index_2, for the additional triangle-to
  //   -tetrahedron mapping.
  //   This function is implemented specifically for the interior
  //   surface between two physical domains, e.g. interior surface
  //   in FSI problem.
  // ----------------------------------------------------------------
  void write_triangle_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<int> &global_node_index,
      const std::vector<int> &global_ele_index_1, 
      const std::vector<int> &global_ele_index_2 );


  // ----------------------------------------------------------------
  // ! write_quadratic_triangle_grid: write the surface mesh described 
  //                                  by quadratic triangle elements
  //                                  with two element index arrays.
  //   The input parameter is identical to the former one, with one
  //   more array, global_ele_index_2, for the additional triangle-to
  //   -tetrahedron mapping.
  //   This function is implemented specifically for the interior
  //   surface between two physical domains, e.g. interior surface
  //   in FSI problem.
  // ----------------------------------------------------------------
  void write_quadratic_triangle_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<int> &global_node_index,
      const std::vector<int> &global_ele_index_1, 
      const std::vector<int> &global_ele_index_2 );


  // ----------------------------------------------------------------
  // ! add_int_PointData : add a point data associated with nodal
  //                       points in the grid_w, which can be either
  //                       a vtkUnstructuredGrid or vtkPolyData.
  //   Input: \para grid_w : the grid vtk object that has been setted
  //                         up with the basic geometry infomation.
  //          \para ptdata : the integer data associated with the grid
  //                         points.
  //          \para dataname : the name of the data to be written.
  // ----------------------------------------------------------------
  void add_int_PointData( vtkPointSet * const &grid_w,
      const std::vector<int> &ptdata, const std::string &dataname );

  void add_double_PointData( vtkPointSet * const &grid_w,
      const std::vector<double> &ptdata, const std::string &dataname );

  void add_Vector3_PointData( vtkPointSet * const &grid_w,
      const std::vector<Vector_3> &ptdata, const std::string &dataname );

  // ----------------------------------------------------------------
  // ! add_int_CellData : add a cell data associated with cells in 
  //                      thegrid_w, which can be either
  //                      a vtkUnstructuredGrid or vtkPolyData.
  //   Input: \para grid_w : the grid vtk object that has been setted
  //                         up with the basic geometry infomation.
  //          \para ptdata : the cell data associated with the grid
  //                         points.
  //          \para dataname : the name of the data to be written.
  // ----------------------------------------------------------------
  void add_int_CellData( vtkPointSet * const &grid_w, 
      const std::vector<int> &cldata, const std::string &dataname );

  void add_double_CellData( vtkPointSet * const &grid_w, 
      const std::vector<double> &cldata, const std::string &dataname );


  // ----------------------------------------------------------------
  // ! write_vtkUnstructuredGrid: write the info in vtkUnstructuredGrid
  //                              to a vtu or vtk file. 
  //   Input: \para filename : filename.{vtu/vtk} is the file to be written
  //          \para grid_w   : vtkUnstructuredGrid object to be written
  //          \para isXML    : flag for vtu (true) or vtk (false)
  // ----------------------------------------------------------------
  void write_vtkPointSet( const std::string &filename, 
      vtkPointSet * const &grid_w, const bool &isXML = true );


  // ================================================================
  // 3. Mesh quality measures
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
  double get_aspect_ratio( const std::vector<double> &coors );


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
  
  // ================================================================
  // 4. TetGen interface
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
  // 5. Tet4 class defines a four node linear tetrahedron object with 
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
      Tet4( const std::vector<double> &nodes );

      Tet4( const std::vector<double> &ctrlPts, 
          const int &ien0, const int &ien1, 
          const int &ien2, const int &ien3 );

      virtual ~Tet4();

      void reset( const std::vector<double> &nodes );

      void reset( const std::vector<double> &ctrlPts,
          const int &ien0, const int &ien1,
          const int &ien2, const int &ien3 );

      // Changes the gindex array only. It is used to determine
      // the face index, e.g. in ElemBC_3D_tet4::resetTriIEN_outwardnormal
      // function.
      void reset( const int &ien0, const int &ien1,
          const int &ien2, const int &ien3 );

      void reset( const std::vector<double> &ctrlPts,
          const IIEN * const &ien_pt, const int &ee );

      double get_aspect_ratio() const;

      double get_volume() const;

      // get the circumscribing sphere's DIAMETER.
      // This is a measure of the mesh size
      double get_diameter() const;

      // Given the face node IEN indices, determine the face id is determined by
      // its opposite node number.
      int get_face_id(const int &n0, const int &n1, const int &n2) const;

      void write_vtu( const std::string &fileName ) const;

      void print_info() const;

    private:
      double pts[12];

      int gindex[4];
  };


  // ================================================================
  // 6. Tetrahedral mesh checker. 
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