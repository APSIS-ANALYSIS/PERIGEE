#include "Map_Node_Index.hpp"

Map_Node_Index::Map_Node_Index( const IGlobal_Part * const &gpart,
    const int &cpu_size, const int &nFunc )
{
  int newnum = 0;
  old_2_new.resize(nFunc);
  new_2_old.resize(nFunc);
  
  std::cout<<"-- generating old2new & new2old index mapping. \n";

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
  std::cout<<std::endl<<"=== Node index mapping generated.\n";
}

Map_Node_Index::Map_Node_Index( const IGlobal_Part * const &gpart,
    const int &cpu_size, const int &n_start, const int &n_end )
{
  int newnum = 0;
  const int nfunc = n_end - n_start;
  old_2_new.resize(nfunc);
  new_2_old.resize(nfunc);

  std::cout<<"\n-- generating old2new & new2old index mapping. \n";

  // old to new mapping
  for(int proc = 0; proc<cpu_size; ++proc)
  {
    for( int nn=0; nn<nfunc; ++nn )
    {
      if( (int) gpart->get_npart(nn+n_start) == proc )
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
  std::cout<<"=== Node index mapping generated.\n";
}

Map_Node_Index::Map_Node_Index( const char * const &fileName )
{
  std::cout<<"-- loading old2new & new2old index mapping from disk. \n";
  std::string fName( fileName );
  fName.append(".h5");

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  old_2_new = h5r -> read_intVector("/", "old_2_new");
  new_2_old = h5r -> read_intVector("/", "new_2_old");

  delete h5r; H5Fclose( file_id );
  
  std::cout<<"-- mapping generated. Memory usage: ";
  SYS_T::print_mem_size( double(old_2_new.size())*2.0*sizeof(int) );
  std::cout<<std::endl<<"=== Node index mapping loaded.\n";
}

Map_Node_Index::~Map_Node_Index()
{
  VEC_T::clean(old_2_new); VEC_T::clean(new_2_old);
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
  
  // file creation
  hid_t file_id = H5Fcreate( fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  h5w -> write_intVector( "old_2_new", old_2_new );
  h5w -> write_intVector( "new_2_old", new_2_old );

  delete h5w; H5Fclose(file_id);
}

// EOF
