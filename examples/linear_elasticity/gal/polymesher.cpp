// ==================================================================
// polymesher.cpp
// This is the mesher I use to generate tetrahedral mesh by directly
// calling tetgen funcitons and the input file is loaded from the 
// .poly format file.
//
// Date: Dec. 23 2016
// ==================================================================

#include "Tet_Tools.hpp"

int main( int argc, char * argv[] )
{
  // ----------------------------------------------------------------
  // Setup the switches for tetgen
  // ----------------------------------------------------------------
  double maxRadRatio = 1.5;
  double minDiheAng = 15.0;
  double maxEdgeSize = 2.0;
  int optimLevel = 2;
  int optimScheme = 7;
  
  char * char_home_dir = getenv("HOME");
  std::string poly_file(char_home_dir);
  poly_file.append("/PERIGEE/input/tet_poly/cube");

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  PetscMPIInt size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  if(size != 1) SYS_T::print_fatal("ERROR: polymesher has to be run in serial! \n");
 
  SYS_T::GetOptionReal("-maxEdgeSize", maxEdgeSize); 
  SYS_T::GetOptionReal("-maxRadRatio", maxRadRatio); 
  SYS_T::GetOptionReal("-minDiheAng", minDiheAng); 
  SYS_T::GetOptionString("-poly_file", poly_file); 
  
  double maxTetVol = std::pow(maxEdgeSize, 3) / (6.0 * std::sqrt(2));

  std::cout<<" =================================================="<<std::endl; 
  std::cout<<" Mesh options: \n";
  std::cout<<"  -- Maximum radius-edge ratio (-maxRadRatio) = "<<maxRadRatio<<std::endl;
  std::cout<<"  -- Minimum dihedral angle (-minDiheAng) = "<<minDiheAng<<" degrees \n";
  std::cout<<"  -- Maximum edge size (-maxEdgeSize) = "<<maxEdgeSize<<'\n';
  std::cout<<"  -- Maximum Tet volume = "<<maxTetVol<<'\n';
  std::cout<<"  -- Optimization level = "<<optimLevel<<'\n';
  std::cout<<"  -- Optimization scheme = "<<optimScheme<<'\n';
  std::cout<<"  -- Input file (-poly_file) = "<<poly_file<<std::endl;

  char switches [300]; 
  int len = sprintf(switches,"pq%.2f/%.1fa%8.3eO%d/%dV",
      maxRadRatio, minDiheAng, maxTetVol, optimLevel, optimScheme);

  SYS_T::print_exit_if( len > 300, "Tetgen option switch is too long. \n");

  std::cout<<"  -- options = "<<switches<<std::endl;
  std::cout<<" =================================================="<<std::endl;

  // ----------------------------------------------------------------
  // Read in the geometry file
  // ----------------------------------------------------------------
  tetgenio in, out;

  std::vector<char> charpolyfile;
  SYS_T::to_char(poly_file, charpolyfile);
  in.load_poly( &charpolyfile[0] );

  // Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
  // do quality mesh generation (q) with a specified quality bound
  // (1.414), and apply a maximum volume constraint (a0.1).
  tetrahedralize( switches, &in, &out);

  // write grid file in .vtu
  std::vector<int> ien;
  const int numcel = out.numberoftetrahedra;
  ien.resize( numcel * 4 );
  for(int ii=0; ii<numcel*4; ++ii) ien[ii] = out.tetrahedronlist[ii] - 1;

  // write output file
  TET_T::tetgenio2vtu( out, "vol" );

  // this output file is for debug only
  TET_T::tetgenio2vtu_windex( out, "vol_index" );

  TET_T::tetgenio2vtp( out, "sur", 1);  
  TET_T::tetgenio2vtp( out, "sur", 2);  
  TET_T::tetgenio2vtp( out, "sur", 3);  
  TET_T::tetgenio2vtp( out, "sur", 4);  
  TET_T::tetgenio2vtp( out, "sur", 5);  
  TET_T::tetgenio2vtp( out, "sur", 6);  

  // Print options again at the bottom of the display
  std::cout<<" =================================================="<<std::endl; 
  std::cout<<" Mesh options: \n";
  std::cout<<"  -- Maximum radius-edge ratio (-maxRadRatio) = "<<maxRadRatio<<std::endl;
  std::cout<<"  -- Minimum dihedral angle (-minDiheAng) = "<<minDiheAng<<" degrees \n";
  std::cout<<"  -- Maximum edge size (-maxEdgeSize) = "<<maxEdgeSize<<'\n';
  std::cout<<"  -- Maximum Tet volume = "<<maxTetVol<<'\n';
  std::cout<<"  -- Optimization level = "<<optimLevel<<'\n';
  std::cout<<"  -- Optimization scheme = "<<optimScheme<<'\n';
  std::cout<<"  -- Input file (-poly_file) = "<<poly_file<<std::endl;
  std::cout<<"  -- options = "<<switches<<std::endl;
  std::cout<<" =================================================="<<std::endl;
  
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
