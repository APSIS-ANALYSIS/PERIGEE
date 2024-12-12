#ifndef VIS_TOOLS_HPP
#define VIS_TOOLS_HPP
// ==================================================================
// Vis_Tools.hpp
//
// This is a set of visualization tools. The original version of these
// tools are written in VTK_Writer class. We separate these tools from
// the VTK_Writer because (1) these tools are logically independent 
// from the data in VTK_Writer; (2) these tools can be reused in other
// VTK_Writers, such as the VTK_Writer_VTK, VTK_Writer_Solids, etc.
//
// This set of tools have a namespace VIS_T.
//
// It is built based on the VTK library.
//
// Date: July 12 2016
// ==================================================================
#include "FEAElement.hpp"
#include "HDF5_Tools.hpp"

#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkQuad.h"
#include "vtkHexahedron.h"
#include "vtkTetra.h"
#include "vtkQuadraticTetra.h"
#include "vtkTriQuadraticHexahedron.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"

namespace VIS_T
{
  // ================================================================
  // The 1st set of tools contains
  //    --- writeGridObject      ===> vtu/vtk
  //    --- writepvtuFile        ===> pvtu
  //    --- writepvdFile_Init    ===> pvd
  //    --- writepvdFile         ===> pvd
  //    --- writeVisFile         ===> write vtu/vtk, pvtu, pvd 
  // ================================================================
  // ----------------------------------------------------------------
  // ! WriteGridObject : Given the girdData that carries the 
  //   visualization infomation, this function writes the vtu/vtk
  //   files on disk. The file name will be 
  //              baseName_pxxxx.vtk/vtu
  //   where xxxx is the rank of the CPU.
  //   If ixXML flag is false, this function will write vtk file.
  // ----------------------------------------------------------------
  void writeGridObject( vtkUnstructuredGrid * gridData,
      const std::string &baseName, const int &rank, const bool &isXML );
  
  // ----------------------------------------------------------------
  // ! writepvtuFile : write the pvtu file that organize all vtu files
  //   with name: baseName.pvtu
  //   This function will check if isXML == true && size > 1.
  //   The pvtu file will be written by the proc 0.
  // ----------------------------------------------------------------
  void writepvtuFile( vtkUnstructuredGrid * gridData,
      const std::string &baseName, const int &rank, const int &size,
      const bool &isXML);

  // ----------------------------------------------------------------
  // ! writepvdFile_Init : write the top lines of pvd file.
  //   If the baseName.pvd file does not exist, this function will
  //   be called to generate a baseName.pvd file with initialization
  //   lines.
  // ----------------------------------------------------------------
  void writepvdFile_Init( const std::string &pvdFName );

  // ----------------------------------------------------------------
  // ! writepvdFile : If baseName.pvd does not exist, the writepvdFile_Init
  //   will be called to generate a pvd file. Then the pvtu/vtu files
  //   will be listed with corresponding time info. 
  //   This function will work if isXML == ture && rank == 0 (size does
  //   not have to be greater than 1).
  // ----------------------------------------------------------------
  void writepvdFile( const std::string &baseName,
      const std::string &pre_pvtuname, const int &rank, const int &size,
      const double &sol_time, const bool &isXML );

  // ----------------------------------------------------------------
  // ! writeVisFile : Given the vtkUnstructuredGrid data is completed,
  //   this function is a driver that will call 
  //                   writeGridObject
  //                   writepvtuFile
  //                   writepvdFile
  //   to write all files.
  //   \para girdData : the prepared vtkGrid data
  //   \para baseName : the base name for the vis output (e.g. SOL_)
  //   \para pre_pvtuname : the pvtu's name (e.g. SOL_9xxxxxxxx)
  //   \para rank  : MPI rank
  //   \para size  : MPI size
  //   \para sol_time : the actual time to be written
  //   \para ixXML : flag for vtk/vtu choice
  // ----------------------------------------------------------------
  void writeVisFile( vtkUnstructuredGrid * gridData, 
      const std::string &baseName, const std::string &pre_pvtuname,
      const int &rank, const int &size, 
      const double &sol_time, const bool &isXML );


  // ================================================================
  // The 2nd set of tools contains 
  //    --- setHexelem      ===> Insert vtkHexahedron to gird
  //    --- setQuadelem     ===> Insert vtkQuad to grid
  // ================================================================
  // ----------------------------------------------------------------
  // ! setHexelem: build Hexahedron element connectivity with the
  //   vtkPoints.  There are 8 points in this element.
  //   \para ptoffset: offset for this element
  //   \gridData: the Grid object that stores the connectivity
  // ----------------------------------------------------------------
  void setHexelem( const int &ptoffset, vtkUnstructuredGrid * gridData );

  // ----------------------------------------------------------------
  // ! setHexelem: build Hexahedron element connectivity using the
  //   specified point index.
  // ----------------------------------------------------------------
  void setHexelem( const int &ptid0, const int &ptid1,
      const int &ptid2, const int &ptid3, const int &ptid4,
      const int &ptid5, const int &ptid6, const int &ptid7,
      vtkUnstructuredGrid * gridData );

  // ----------------------------------------------------------------
  // ! setTriQuadHexelem: build triquadratic Hex element connectivity
  //   with the vtkPoints.  There are 27 points in this element.
  //   \para ptoffset: offset for this element
  //   \gridData: the Grid object that stores the connectivity
  // ----------------------------------------------------------------
  void setTriQuadHexelem( const int &ptoffset, vtkUnstructuredGrid * gridData );

  // ----------------------------------------------------------------
  // ! setTriQuadHexelem: build triquadratic Hex element connectivity
  //   using the specified point index.
  // ----------------------------------------------------------------
  void setTriQuadHexelem( const int &ptid0, const int &ptid1,
      const int &ptid2, const int &ptid3, const int &ptid4,
      const int &ptid5, const int &ptid6, const int &ptid7,
      const int &ptid8, const int &ptid9, const int &ptid10,
      const int &ptid11, const int &ptid12, const int &ptid13,
      const int &ptid14, const int &ptid15, const int &ptid16,
      const int &ptid17, const int &ptid18, const int &ptid19,
      const int &ptid20, const int &ptid21, const int &ptid22,
      const int &ptid23, const int &ptid24, const int &ptid25,
      const int &ptid26, vtkUnstructuredGrid * gridData );
 
  // ----------------------------------------------------------------
  // ! setTetraelem: build Tetra element connectivity with vtkPoints
  //   in element. There are 4 points in this element.
  //   \para ptoffset: the offset for this element
  //   \para gridData: the Grid object that stores the connectivity
  // ----------------------------------------------------------------
  void setTetraelem( const int &ptoffset, vtkUnstructuredGrid * gridData );

  // ----------------------------------------------------------------
  // ! setTetraelem: build Tetra element connectivity using the specified
  //   point index.
  // ----------------------------------------------------------------
  void setTetraelem( const int &ptid0, const int &ptid1, const int &ptid2,
      const int &ptid3, vtkUnstructuredGrid * gridData );

  // ----------------------------------------------------------------
  // ! setQuadTetraelem: build quadratic Tetra element connectivity
  //   with the vtkPoints. There are 10 points in the element.
  //   \para ptoffset: the offset for this element
  //   \para gridData: the Grid object that stores the connectivity
  // ----------------------------------------------------------------
  void setQuadTetraelem( const int &ptoffset, vtkUnstructuredGrid * gridData );

  // ----------------------------------------------------------------
  // ! setQuadTetraelem: build quadratic Tetra element connectivity 
  //   using the specified point index.
  // ----------------------------------------------------------------
  void setQuadTetraelem( const int &ptid0, const int &ptid1,
      const int &ptid2, const int &ptid3, const int &ptid4,
      const int &ptid5, const int &ptid6, const int &ptid7,
      const int &ptid8, const int &ptid9,
      vtkUnstructuredGrid * gridData );

  // ================================================================
  // The 3rd set of tools contain
  //    --- read_epart     ===> Read epart.h5
  // ================================================================
  // --------------------------------------------------------------
  // ! read_epart: read the element partition file, usually named as
  //               "epart.h5"
  // \para epart_file: epart.h5, i.e., the element partition from 
  //                   analysis code's preprocessor
  // \para esize: the total element, should be read from global info
  // \output elem_part: the int array with size esize, is the "part" in
  //                    "epart.h5". Users are responsible for allocating
  //                    the array with size esize, and deleting it after
  //                    usage.
  // --------------------------------------------------------------
  void read_epart( const std::string &epart_file, const int &esize,
      std::vector<int> &elem_part );
}

#endif
