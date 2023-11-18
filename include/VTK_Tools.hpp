#ifndef VTK_TOOLS_HPP
#define VTK_TOOLS_HPP
// ============================================================================
// VTK_Tools.hpp
//
// This is a suite of function tools that assist some fundamental operations for
// VTK files.
// 
// Date Created: Aug. 12 2023
// ============================================================================
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"

#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLGenericDataObjectReader.h"

namespace VTK_T
{
  // ================================================================
  // ===> 1. The first set of tools assists READING volumetric mesh 
  //         from .vtu file and surface mesh from .vtp file.  
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
  // ! read_grid: read a generic mesh from either a .vtp or a .vtu file.
  //              The output is identical to the read_vtp_grid and read_vtu_grid
  //              functions.
  //   Note: it will return 1 if the file is of vtp type,
  //                 return 2 if the file is of vtu type,
  //                 return 0 if the file is unindentified.
  // ----------------------------------------------------------------
  int read_grid( const std::string &filename,
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

  std::vector<Vector_3> read_Vector3_PointData( const std::string &filename,
      const std::string &dataname );

  // ----------------------------------------------------------------
  // ! read_num_pt: read the number of points from either a .vtp or a .vtu file.
  // Input: \para filename the vtk file name
  // Output: the number of points in the file
  // ----------------------------------------------------------------
  int read_num_pt( const std::string &filename );

  // ----------------------------------------------------------------
  // ! read_num_cl: read the number of cells from either a .vtp or a .vtu file.
  // Input: \para filename the vtk file name
  // Output: the number of cells in the file
  // ----------------------------------------------------------------
  int read_num_cl( const std::string &filename );

  // ================================================================
  // ===> 2. The second set of tools assists WRITING volumetric mesh 
  //         to .vtu file and surface mesh to .vtp file.  
  // ================================================================
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
  // ! write_vtkPointSet : write the info in vtkUnstructuredGrid or
  //                       vtkXMLPolyData to a vtu or vtk file. 
  //   Input: \para filename : filename.{vtu/vtk} is the file to be written.
  //          \para grid_w   : vtkUnstructuredGrid or vtkXMLPolyData object 
  //                           to be written.
  //          \para isXML    : flag for vtu (true) or vtk (false).
  // ----------------------------------------------------------------
  void write_vtkPointSet( const std::string &filename, 
      vtkPointSet * const &grid_w, const bool &isXML = true );
}

#endif
