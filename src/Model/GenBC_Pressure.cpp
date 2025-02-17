#include "GenBC_Pressure.hpp"

GenBC_Pressure::GenBC_Pressure( const std::string &lpn_filename, const double &init_time )
{
  SYS_T::file_check( lpn_filename ); // make sure the file is on the disk

  std::ifstream reader;
  reader.open( lpn_filename.c_str(), std::ifstream::in );

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

  if( bc_type.compare("Pressure") == 0 || bc_type.compare("PRESSURE") == 0 )
  {
    coef_a.resize(num_ebc); coef_b.resize(num_ebc);
    num_of_mode.resize(num_ebc); w.resize(num_ebc); period.resize(num_ebc);
  }
  else
    SYS_T::print_fatal( "GenBC_Pressure Error: BC type in %s should be Pressure.\n", lpn_filename.c_str() );

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

        SYS_T::print_fatal_if( face_id != ebc_id, "GenBC_Pressure Error: ebc in %s should be listed in ascending order.\n", lpn_filename.c_str() );

        sstrm >> num_of_mode[ebc_id];
        sstrm >> w[ebc_id];
        sstrm >> period[ebc_id];

        sstrm.clear();
        break;
      }
    }

    // Check the compatibility of period and w. If the difference
    // is larger than 0.01, print a warning message
    if( std::abs(2.0 * MATH_T::PI / period[ebc_id] - w[ebc_id] ) >= 0.01 ) SYS_T::commPrint( "\nGenBC_Pressure WARNING: ebc_id %d incompatible period and w, \n2xpi/period = %e and w = %e.\n", ebc_id, 2.0*MATH_T::PI/period[ebc_id], w[ebc_id] );

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

    coef_a[ebc_id].shrink_to_fit();

    if( static_cast<int>(coef_a[ebc_id].size()) != num_of_mode[ebc_id]+1 ) SYS_T::print_fatal( "GenBC_Pressure Error: ebc_id %d a-coefficients in %s incompatible with the given number of modes.\n", ebc_id, lpn_filename.c_str() );

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

    coef_b[ebc_id].shrink_to_fit();

    if( static_cast<int>(coef_b[ebc_id].size()) != num_of_mode[ebc_id]+1 ) SYS_T::print_fatal( "GenBC_Pressure Error: ebc_id %d b-coefficients in %s incompatible with the given number of modes.\n", ebc_id, lpn_filename.c_str() );
  }

  // Finish reading the file and close it
  reader.close();

  // Initialize P0 by setting it to be get_P at time = 0.0
  P0.resize( num_ebc );
  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
    P0[ebc_id] = get_P( ebc_id, 0.0, 0.0, init_time );
}

void GenBC_Pressure::print_info() const
{
  SYS_T::commPrint("===> GenBC_Pressure: \n");

  for(int ii=0; ii<num_ebc; ++ii)
  {
    SYS_T::commPrint("     ebc_id = %d", ii);
    SYS_T::commPrint(" w = %e, period =%e \n", w[ii], period[ii]);
    SYS_T::commPrint("     a[0] + Sum{ a[i] cos(i x w x t) + b[i] sin(i x w x t) }, for i = 1,...,%d. \n", num_of_mode[ii]);
    for(int jj=0; jj<=num_of_mode[ii]; ++jj)
      SYS_T::commPrint("     i = %d, a = %e, b = %e \n", jj, coef_a[ii][jj], coef_b[ii][jj]);
  }
}

double GenBC_Pressure::get_P( const int &ebc_id, const double &dot_Q, const double &Q,
       const double &time ) const
{
  double PP = coef_a[ebc_id][0];
  for( int ii = 1; ii <= num_of_mode[ebc_id]; ++ii )
    PP += coef_a[ebc_id][ii] * cos( ii*w[ebc_id]*time ) + coef_b[ebc_id][ii] * sin( ii*w[ebc_id]*time );

  return PP;
}

// EOF
