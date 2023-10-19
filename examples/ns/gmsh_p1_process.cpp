// ==================================================================
// gmsh_p1_process.cpp
//
// Code that handles the gmsh file generated using linear tet elem.
//
// Date Created: July 1 2017
// Author: Ju Liu
// ==================================================================
#include "Gmsh_FileIO.hpp"

int main( int argc, char * argv[] )
{
  std::string gmshFile = "cylinder.msh";
  int num_outlet = 1;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  SYS_T::GetOptionString("-gmsh_file", gmshFile);
  SYS_T::GetOptionInt("-num_outlet", num_outlet);
  
  std::cout<<" -gmsh_file: "<<gmshFile<<std::endl;
  std::cout<<" -num_outlet: "<<num_outlet<<std::endl;

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  GIO -> print_info();

  // Gmsh file should contain only one 3D domain
  const int vol_id = 0;

  GIO -> write_vtp("wall_vol", 0, vol_id, true); // assumed to be wall vtp
  GIO -> write_vtp("inflow_vol_000", 1, vol_id, true); // assumed to be inlet vtp
  
  for(int ii=0; ii<num_outlet; ++ii)
  {
    std::string filename {"outflow_vol_"};
    std::string idx = std::to_string(ii);
    if(ii < 10)
    {
      filename += "00";
      filename += idx;
    }
    else if(ii > 9 && ii < 100)
    {
      filename += "0";
      filename += idx;
    }
    else if(ii > 99 && ii < 1000)
      filename += idx;
    else
      SYS_T::print_fatal("Too many outlets.");
    
    GIO -> write_vtp(filename, 2+ii, vol_id, true);
  }

  const std::string wmname("whole_vol");
  const bool isXML = true;
  GIO -> write_vtu( wmname, isXML );

  delete GIO; 
  PetscFinalize();
  return 0;
}

// EOF
