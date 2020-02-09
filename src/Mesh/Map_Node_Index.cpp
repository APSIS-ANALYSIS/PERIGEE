#include "Map_Node_Index.hpp"

Map_Node_Index::Map_Node_Index( const class IGlobal_Part * const &gpart,
    const int &cpu_size, const int nFunc )
{
  int newnum = 0;
  old_2_new.resize(nFunc);
  new_2_old.resize(nFunc);
  
  std::cout<<"\n-- generating old2new & new2old index mapping. \n";

  // old to new mapping
  for(int proc = 0; proc<cpu_size; ++proc)
  {
    for( int nn=0; nn<nFunc; ++nn )
    {
      if( (int) gpart->get_npart(nn) == proc )
      {
        old_2_new[nn] = newnum;
        new_2_old[newnum] = nn;
        newnum += 1;
      }
    }
  }
  VEC_T::shrink2fit(old_2_new);
  VEC_T::shrink2fit(new_2_old);
  
  std::cout<<"-- mapping generated. Memory usage: ";
  SYS_T::print_mem_size( double(old_2_new.size())*2.0*sizeof(int) );
  std::cout<<std::endl;
  std::cout<<"=== Node index mapping generated. \n \n";
}

Map_Node_Index::~Map_Node_Index()
{
  VEC_T::clean(old_2_new); VEC_T::clean(new_2_old);
  std::cout<<"-- Map_Node_Index deleted. \n";
}

void Map_Node_Index::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"========================== "<<std::endl;
  std::cout<<" ii - old_2_new - new_2_old "<<std::endl;
  for(int ii=0; ii<(int) old_2_new.size(); ++ii)
    std::cout<<ii<<'\t'<<get_old2new(ii)<<'\t'<<get_new2old(ii)<<std::endl;
  std::cout<<std::endl;
  std::cout<<" Memory usage: "<<old_2_new.size() * 2 * sizeof(int)<<" bytes."<<std::endl;
  std::cout<<"========================== "<<std::endl;
}

void Map_Node_Index::write_hdf5( const char * const &fileName ) const
{
  std::string fName(fileName);
  fName.append(".h5");
  
  hid_t file_id, dataspace_id;
  hsize_t dim[1];

  // file creation
  file_id = H5Fcreate( fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  // dataspace
  dim[0] = old_2_new.size();
  dataspace_id = H5Screate_simple(1, dim, NULL);

  // dataset
  hid_t setid1;
  setid1 = H5Dcreate(file_id, "old_2_new", H5T_NATIVE_INT, dataspace_id,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(setid1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &old_2_new[0]);
  H5Dclose(setid1);

  hid_t setid2;
  setid2 = H5Dcreate(file_id, "new_2_old", H5T_NATIVE_INT, dataspace_id,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(setid2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &new_2_old[0]);
  H5Dclose(setid2);

  H5Sclose(dataspace_id);
  H5Fclose(file_id);
}


// EOF
