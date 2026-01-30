#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>

#include "Sys_Tools.hpp"
#include "HDF5_Tools.hpp"
#include "Vis_Tools.hpp"

// This ReadNodeMapping function is from PostVectSolution class
void ReadNodeMapping_1( const std::string &node_mapping_file,
    const char * const &mapping_type, const int &node_size,
    int * const &nodemap )
{
  const std::vector<int> temp_nodemap = HDF5_T::read_intVector( node_mapping_file.c_str(), "/", mapping_type );

  SYS_T::print_fatal_if(int( temp_nodemap.size() ) != node_size, "Error: PostVectSolution the allocated array has wrong size! \n");

  for(int ii=0; ii<node_size; ++ii) nodemap[ii] = temp_nodemap[ii];
}

// This ReadNodeMapping function is from vis_wss_*.cpp
std::vector<int> ReadNodeMapping_2( const char * const &node_mapping_file,
    const char * const &mapping_type, const int &node_size )
{
  hid_t file_id = H5Fopen(node_mapping_file, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t data_id = H5Dopen(file_id, mapping_type, H5P_DEFAULT);

  hid_t data_space = H5Dget_space( data_id );
  hid_t data_rank = H5Sget_simple_extent_ndims( data_space );

  if( data_rank != 1)
  {
    SYS_T::commPrint("Error: the node mapping file has wrong format. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  hsize_t * data_dims = new hsize_t [1];

  H5Sget_simple_extent_dims( data_space, data_dims, NULL );

  hid_t mem_space = H5Screate_simple(data_rank, data_dims, NULL);

  hsize_t dSize = data_dims[0];

  if( int(dSize) != node_size )
  {
    SYS_T::commPrint("Error: the allocated array has wrong size! \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  std::vector<int> out(node_size, -1);

  H5Dread( data_id, H5T_NATIVE_INT, mem_space, data_space, H5P_DEFAULT, &out[0] );

  delete [] data_dims;
  H5Sclose( mem_space );
  H5Sclose(data_space);
  H5Dclose(data_id);
  H5Fclose(file_id);

  return out;
}

int main()
{
  const std::string node_mapping_file = "node_mapping.h5";
  
  const int nNode = 3456;

  std::vector<int> new2old_map_1(nNode);
  ReadNodeMapping_1( node_mapping_file, "new_2_old", nNode, new2old_map_1.data() );

  std::vector<int> new2old_map_2 = ReadNodeMapping_2( node_mapping_file.c_str(), "new_2_old", nNode );

  std::cout<<"Map_1: size = " << new2old_map_1.size() << "\n";
  VEC_T::print(new2old_map_1);
  std::cout<<"\n\nMap_2: size = " << new2old_map_2.size() << "\n";
  VEC_T::print(new2old_map_2);

  if(VEC_T::is_equal(new2old_map_1, new2old_map_2, 0)) std::cout<<"good!\n";
  else std::cout<<"bad\n";


  std::vector<int> old2new_map_1(nNode);
  ReadNodeMapping_1( node_mapping_file, "old_2_new", nNode, old2new_map_1.data() );

  std::vector<int> old2new_map_2 = ReadNodeMapping_2( node_mapping_file.c_str(), "old_2_new", nNode );

  std::cout<<"Map_1: size = " << old2new_map_1.size() << "\n";
  VEC_T::print(old2new_map_1);
  std::cout<<"\n\nMap_2: size = " << old2new_map_2.size() << "\n";
  VEC_T::print(old2new_map_2);

  if(VEC_T::is_equal(old2new_map_1, old2new_map_2, 0)) std::cout<<"good!\n";
  else std::cout<<"bad\n";

  std::vector<int> old2new_map_3 = VIS_T::readNodeMapping( node_mapping_file, "old_2_new", nNode );

  if(VEC_T::is_equal(old2new_map_1, old2new_map_3, 0)) std::cout<<"good!\n";
  else std::cout<<"bad\n";

  return EXIT_SUCCESS;
}

// EOF
