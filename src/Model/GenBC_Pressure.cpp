#include "GenBC_Pressure.hpp"

GenBC_Pressure::GenBC_Pressure( const char * const &pres_file_name )
{
  std::ifstream reader;
  reader.open( pres_file_name, std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;
  std::string bc_type;

  // The first non-commented, non-empty line should be
  // Pressure num_ebc
  while( std::getline(reader, sline) )
  {
    if( sline[0] !='#' && !sline.empty() )
    {
      sstrm.str(sline);
      sstrm >> bc_type;
      sstrm >> num_ebc;
      sstrm.clear();
      break;
    }
  }

  if( num_ebc != part_ebc->get_num_ebc() )
    SYS_T::print_fatal( "PLocAssem_Tet_CMM_GenAlpha Error: num_ebc in %s is incompatible with ALocal_EBC", "pressure_fourier_series.txt" );

  if( bc_type.compare("Pressure") == 0 || bc_type.compare("PRESSURE") == 0 )
  {
    coef_a.resize(num_ebc); coef_b.resize(num_ebc);
    num_of_mode.resize(num_ebc); w.resize(num_ebc); period.resize(num_ebc);
  }
  else
    SYS_T::print_fatal( "PLocAssem_Tet_CMM_GenAlpha Error: inlet BC type in %s should be Pressure.\n", "pressure_fourier_series.txt" );

  // Read in num_of_mode, w, period, coef_a, and coef_b per ebc
  for(int ebc_id=0; ebc_id<num_ebc; ++ebc_id)
  {
    while( std::getline(reader, sline) )
    {
      // face_id num_of_mode w period
      if( sline[0] !='#' && !sline.empty() )
      {
        sstrm.str(sline);
        int face_id;
        sstrm >> face_id;

        if( face_id != ebc_id )
          SYS_T::print_fatal( "PLocAssem_Tet_CMM_GenAlpha Error: ebc in %s should be listed in ascending order.\n", "pressure_fourier_series.txt" );

        sstrm >> num_of_mode[ebc_id];
        sstrm >> w[ebc_id];
        sstrm >> period[ebc_id];

        sstrm.clear();
        break;
      }
    }

    // Check the compatibility of period and w. If the difference
    // is larger than 0.01, print a warning message
    if( std::abs(2.0 * MATH_T::PI / period[ebc_id] - w[ebc_id] ) >= 0.01 )
      SYS_T::commPrint( "\nPLocAssem_Tet_CMM_GenAlpha WARNING: ebc_id %d incompatible period and w, \n2xpi/period = %e and w = %e.\n", ebc_id, 2.0*MATH_T::PI/period[ebc_id], w[ebc_id] );

    coef_a[ebc_id].clear(); coef_b[ebc_id].clear();

    while( std::getline(reader, sline) )
    {
      // coef_a
      if( sline[0] !='#' && !sline.empty() )
      {
        sstrm.str(sline);
        double temp_coef;
        while( sstrm >> temp_coef ) coef_a[ebc_id].push_back( temp_coef );

        sstrm.clear();
        break;
      }
    }

    VEC_T::shrink2fit( coef_a[ebc_id] );

    if( static_cast<int>(coef_a[ebc_id].size()) != num_of_mode[ebc_id]+1 )
      SYS_T::print_fatal( "PLocAssem_Tet_CMM_GenAlpha Error: ebc_id %d a-coefficients in %s incompatible with the given number of modes.\n", ebc_id, "pressure_fourier_series.txt" );

    while( std::getline(reader, sline) )
    {
      // coef_b
      if( sline[0] !='#' && !sline.empty() )
      {
        sstrm.str(sline);
        double temp_coef;
        while( sstrm >> temp_coef ) coef_b[ebc_id].push_back( temp_coef );

        sstrm.clear();
        break;
      }
    }

    VEC_T::shrink2fit( coef_b[ebc_id] );

    if( static_cast<int>(coef_b[ebc_id].size()) != num_of_mode[ebc_id]+1 )
      SYS_T::print_fatal( "PLocAssem_Tet_CMM_GenAlpha Error: ebc_id %d b-coefficients in %s incompatible with the given number of modes.\n", ebc_id, "pressure_fourier_series.txt" );
  }

  // Finish reading the file and close it
  reader.close();
}

GenBC_Pressure::~GenBC_Pressure()
{}

void GenBC_Pressure::print_info() const
{}

double GenBC_Pressure::get_P( const int &ii, const double &dot_Q, const double &Q,
       const double &time ) const
{}

double GenBC_Pressure::get_P0( const int &ii ) const
{}

void GenBC_Pressure::reset_initial_sol( const int &ii, const double &in_Q_0,
        const double &in_P_0, const double &curr_time, const bool &is_restart )
{}

// EOF
