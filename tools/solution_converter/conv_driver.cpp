// ==================================================================
// conv_driver.cpp
//
// This is a converter for the PDNSolution file. If the same mesh
// is partitioned in two different ways, this driver can help convert
// the solution vector numbering to make the solution compatible with
// the new partitioning.
//
// A common usage of this file is to convert the solution for problems
// run with different part files from different -cpu_size argument.
//
// Example:
// ./sol_converter -old_nmap ../build_ns/node_mapping.h5 
//                 -new_nmap ../build_ns2/node_mapping.h5 
//                 -sol_name SOL_900000001 
//                 -out_name NEW_900000001
//
//  SOL_900000001 is compatible with old nmap nodes;
//  NEW_900000001 is compatible with new nmap nodes.
//
// Author: Ju Liu
// Date: Mar. 13 2019
// ==================================================================
#include "Sys_Tools.hpp"
#include "HDF5_Reader.hpp"

int main( int argc, char * argv[] )
{
  std::string old_nmap("old_node_mapping.h5");
  std::string new_nmap("new_node_mapping.h5");
  std::string sol_name("SOL_900000000");
  std::string out_name("NEW_900000000");
  

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const PetscMPIInt size = SYS_T::get_MPI_size();
  
  SYS_T::print_fatal_if(size!=1,"ERROR: converter is a serial routine! \n");

  SYS_T::GetOptionString("-old_nmap", old_nmap);
  SYS_T::GetOptionString("-new_nmap", new_nmap);
  SYS_T::GetOptionString("-sol_name", sol_name);
  SYS_T::GetOptionString("-out_name", out_name);
  
  std::cout<<"==== /Command Line Arguments ===="<<std::endl;
  std::cout<<" -old_nmap: "<<old_nmap<<std::endl;
  std::cout<<" -new_nmap: "<<new_nmap<<std::endl;
  std::cout<<" -sol_name: "<<sol_name<<std::endl;
  std::cout<<" -out_name: "<<out_name<<std::endl;

  SYS_T::file_check(old_nmap);
  SYS_T::file_check(new_nmap);
  SYS_T::file_check(sol_name);

  // Obtain the new_to_old vector from the old mapping
  hid_t h5id_onm = H5Fopen(old_nmap.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * onm_h5r = new HDF5_Reader( h5id_onm );

  std::vector<int> old_map = onm_h5r-> read_intVector( "/", "old_2_new" );

  delete onm_h5r; H5Fclose(h5id_onm);

  // Obtain the old_to_new vector from the new mapping
  hid_t h5id_nnm = H5Fopen(new_nmap.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * nnm_h5r = new HDF5_Reader( h5id_nnm );

  std::vector<int> new_map = nnm_h5r-> read_intVector( "/", "old_2_new" );

  delete nnm_h5r; H5Fclose(h5id_nnm);

  // Check the length of the two vectors, if they match, we assign
  // the lenght to the value of nFunc
  SYS_T::print_fatal_if( old_map.size() != new_map.size(),
     "Error: The two node mapping files node index lengths do not match.\n" );

  const int nFunc = static_cast<int>( old_map.size() );

  // Read the PETSc solution
  Vec sol_temp;
  VecCreate(PETSC_COMM_SELF, &sol_temp);
  VecSetType(sol_temp, VECSEQ);

  PetscViewer viewer_read;
  PetscViewerBinaryOpen(PETSC_COMM_SELF, sol_name.c_str(), FILE_MODE_READ, &viewer_read);
  VecLoad(sol_temp, viewer_read);
  PetscViewerDestroy(&viewer_read);
  
  // Check the solution length and determine dof
  PetscInt sol_size;
  VecGetSize(sol_temp, &sol_size);
  SYS_T::print_fatal_if( sol_size % nFunc != 0, 
      "Error: Solution file is incompatible with the node mapping file. \n" );

  const int dof = sol_size / nFunc;

  std::cout<<"Solution node number is "<<nFunc<<" with "<<dof<<" dof(s).\n";

  // Convert sol into array
  double * sol_array = new double [sol_size];
  double * array_temp;
  VecGetArray(sol_temp, &array_temp);

  for(int ii=0; ii<sol_size; ++ii) sol_array[ii] = array_temp[ii];

  VecRestoreArray(sol_temp, &array_temp);
  VecDestroy(&sol_temp);

  // Converting index
  Vec sol_new;
  VecCreate(PETSC_COMM_SELF, &sol_new);
  VecSetSizes(sol_new, PETSC_DECIDE, sol_size);
  VecSetType(sol_new, VECSEQ);
  
  for(int ii=0; ii<nFunc; ++ii)
  {
    const int idx_old = old_map[ii];
    const int idx_new = new_map[ii];
    for(int jj=0; jj<dof; ++jj) VecSetValue(sol_new, idx_new*dof+jj, sol_array[idx_old*dof+jj], INSERT_VALUES );
  }

  VecAssemblyBegin(sol_new); VecAssemblyEnd(sol_new);

  // Write new solution to disk
  PetscViewer viewer;
  PetscViewerCreate(PETSC_COMM_SELF, &viewer);
  PetscViewerSetType(viewer, PETSCVIEWERBINARY);
  PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
  PetscViewerBinarySkipInfo(viewer);
  PetscViewerFileSetName(viewer, out_name.c_str());
  VecView(sol_new, viewer);
  PetscViewerDestroy(&viewer);

  // new_array is no more needed
  VecDestroy(&sol_new);
  delete [] sol_array; sol_array = nullptr;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
