// ==================================================================
// gmsh_process.cpp
//
// Code that handles the gmsh file and write data into HDF5 files.
//
// Date Created: Dec. 6 2017
// Author: Ju Liu
// ==================================================================
#include "Gmsh_FileIO.hpp"

int main( int argc, char * argv[] )
{
  std::string gmshFile = "tissue.msh";

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  SYS_T::GetOptionString("-gmsh_file", gmshFile);
  std::cout<<" -gmsh_file: "<<gmshFile<<std::endl;

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  GIO -> print_info();

  const int vol_idx = 0;

  std::vector<int> face_idx; face_idx.clear();
  face_idx.push_back( 0 );
  face_idx.push_back( 1 );
  face_idx.push_back( 2 );
  face_idx.push_back( 3 );
  face_idx.push_back( 4 );
  face_idx.push_back( 5 );

  std::vector<int> needmap; needmap.clear();
  needmap.push_back( 0 ); // Only need top for ebc
  needmap.push_back( 1 );
//  needmap.push_back( 2 );
//  needmap.push_back( 3 );
//  needmap.push_back( 4 );
//  needmap.push_back( 5 );

  GIO->write_tet_h5( vol_idx, face_idx, needmap );

  delete GIO; 
  PetscFinalize();
  return 0;
}

// EOF
