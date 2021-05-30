// ============================================================================
// wall_solver.cpp
//
// Triangle element based finite element code for wall mechanics.
//
// ============================================================================

int main( int argc, char *argv[] )
{
  // We assume that a 3D solver has been called (to generate the wall traction)
  // and a suite of command line arguments has been saved to disk
  hid_t solver_cmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  const double wall_density = cmd_h5r -> read_doubleScalar("/", "wall_density");
  const double wall_poisson = cmd_h5r -> read_doubleScalar("/", "wall_poisson");
  const double wall_kappa   = cmd_h5r -> read_doubleScalar("/", "wall_kappa");

  delete cmd_h5r; H5Fclose(prepcmd_file);


  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  PetscFinalize();
  return EXIT_SUCCESS;
}







// EOF
