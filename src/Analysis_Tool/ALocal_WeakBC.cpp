#include "ALocal_WeakBC.hpp"

ALocal_WeakBC::ALocal_WeakBC( const std::string &fileBaseName, const int &cpu_rank)
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/weak");

  weakbc_type = h5r -> read_intScalar( gname.c_str(), "weak_bc_type" );

  if (weakbc_type > 0)
  {
    num_LD = h5r -> read_intScalar( gname.c_str(), "num_LD_weak" ); // dependent from NBC

    num_local_cell = h5r -> read_intScalar( gname.c_str(), "num_local_cell_weak" ); // dependent from EBC

    cell_nLocBas = h5r -> read_intScalar( gname.c_str(), "cell_nLocBas" );  // shared with EBC

    local_cell_ien = h5r -> read_intVector( gname.c_str(), "local_cell_ien_weak" ); // dependent from EBC

    LDN = h5r -> read_intVector( gname.c_str(), "LDN_weak" ); // dependent from NBC

    LD_xyz = h5r -> read_doubleVector( gname.c_str(), "LD_xyz");  // dependent from EBC
  }

  if(weakbc_type == 2)
  {
    const std::vector<double> Q_vec = h5r -> read_doubleVector( gname.c_str(), "wall_rotation_matrix" );  // dependent from ringBC

    Q.resize( num_LD );
    for(int ii{0}; ii < num_LD; ++ii)
    {
      Q[ii] = Tensor2_3D( Q_vec[9*ii+0], Q_vec[9*ii+1], Q_vec[9*ii+2],
                          Q_vec[9*ii+3], Q_vec[9*ii+4], Q_vec[9*ii+5],
                          Q_vec[9*ii+6], Q_vec[9*ii+7], Q_vec[9*ii+8]  );
    }
  }

}