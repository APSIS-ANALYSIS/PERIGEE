#include "Tissue_prestress.hpp"

Tissue_prestress::Tissue_prestress( 
    const ALocal_Elem * const &locelem, const int &in_nqp_tet, 
    const int &in_cpu_rank, const bool &load_from_file,
    const int &num_layer, const std::string &in_ps_fName )
: cpu_rank( in_cpu_rank ), nlocalele(locelem->get_nlocalele()), nqp(in_nqp_tet),
  ps_fileBaseName( in_ps_fName )
{
  qua_prestress.resize( nlocalele );

  int counter_elem_s = 0;
  for(int ee=0; ee<nlocalele; ++ee)
  {
    if( locelem->get_elem_tag(ee) == 0 ) qua_prestress[ee].clear();
    else if( > 0 && locelem->get_elem_tag(ee) <= num_layer )
    {
      qua_prestress[ee].resize(nqp * 6);
      counter_elem_s += 1;
    }
    else
      SYS_T::print_fatal("Error: element tag is wrong.\n");
  }

  // Assign the value of qua_prestress data structure
  if( counter_elem_s > 0 )
  {
    // Assign the qua_prestress to be all zero.
    std::vector<double> qua_ps_array( counter_elem_s * nqp * 6, 0.0 );

    if( load_from_file == true )
    {
      // If the prestress data exist on disk, read the prestress data;
      const std::string ps_fName = SYS_T::gen_partfile_name( ps_fileBaseName, cpu_rank );
      if( SYS_T::file_exist(ps_fName) )
      {
        hid_t ps_file_id = H5Fopen(ps_fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

        HDF5_Reader * ps_h5r = new HDF5_Reader( ps_file_id );

        const int ps_size = ps_h5r -> read_intScalar("/", "ps_array_size"); 

        SYS_T::print_fatal_if(ps_size != counter_elem_s * nqp * 6, "Error: Tissue_prestress the HDF5 file for prestress is incompatible with the local solid element number.\n");

        qua_ps_array = ps_h5r -> read_doubleVector("/", "prestress");

        delete ps_h5r; H5Fclose(ps_file_id);
      }
      else
        SYS_T::print_fatal("Error: prestress file %s cannot be found.\n", ps_fName.c_str());

      // Load the qua_ps_array to the proper location of qua_prestress 
      int offset = 0;
      for(int ee=0; ee<nlocalele; ++ee)
      {
        if( locelem -> get_elem_tag(ee) > 0 )
        {
          for(int ii=0; ii<6*nqp; ++ii)
            qua_prestress[ee][ii] = qua_ps_array[ offset + ii ];

          offset += nqp * 6;
        }
      }
    } // Finish loading prestress from disk file
  }// Finish initializing the prestress data structure
}

Tissue_prestress::~Tissue_prestress()
{
  for(int ee=0; ee<nlocalele; ++ee) VEC_T::clean( qua_prestress[ee] );

  VEC_T::clean( qua_prestress );
}

std::vector<double> Tissue_prestress::get_prestress( const int &ee ) const
{
  return qua_prestress[ee];
}

std::array<double,6> Tissue_prestress::get_prestress( const int &ee, const int &qua ) const
{
  return {{ qua_prestress[ee][6*qua], qua_prestress[ee][6*qua+1], qua_prestress[ee][6*qua+2],
    qua_prestress[ee][6*qua+3], qua_prestress[ee][6*qua+4], qua_prestress[ee][6*qua+5] }};
}

void Tissue_prestress::add_prestress( const int &ee, const double * const &in_psval )
{
  for(int ii=0; ii<6*nqp; ++ii) qua_prestress[ee][ii] += in_psval[ii];
}

void Tissue_prestress::add_prestress( const int &ee, const int &qua, 
    const double * const &in_psval )
{
  qua_prestress[ee][qua*6  ] += in_psval[0];
  qua_prestress[ee][qua*6+1] += in_psval[1];
  qua_prestress[ee][qua*6+2] += in_psval[2];
  qua_prestress[ee][qua*6+3] += in_psval[3];
  qua_prestress[ee][qua*6+4] += in_psval[4];
  qua_prestress[ee][qua*6+5] += in_psval[5];
}

void Tissue_prestress::print_info() const
{
  for(int ee=0; ee<nlocalele; ++ee)
  {
    std::cout<<"ee: "<<ee<<'\t';
    VEC_T::print( qua_prestress[ee] );
    std::cout<<std::endl;
  }
}

void Tissue_prestress::write_prestress_hdf5() const
{
  // Prepare the data to be reccorded
  std::vector<double> qua_prestress_array;
  qua_prestress_array.clear();

  for(int ee=0; ee<nlocalele; ++ee)
  {
    if( qua_prestress[ee].size() > 0 )
      VEC_T::insert_end(qua_prestress_array, qua_prestress[ee]);
  }

  // Record to h5 file
  const std::string fName = SYS_T::gen_partfile_name( ps_fileBaseName, cpu_rank );

  hid_t file_id = H5Fcreate(fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  h5w -> write_intScalar( "ps_array_size", qua_prestress_array.size() );

  if( qua_prestress_array.size() > 0 )
    h5w -> write_doubleVector( file_id, "prestress", qua_prestress_array );

  delete h5w; H5Fclose( file_id );
}

// EOF
