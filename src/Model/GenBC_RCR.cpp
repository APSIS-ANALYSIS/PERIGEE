#include "GenBC_RCR.hpp"

GenBC_RCR::GenBC_RCR( const std::string &lpn_filename, const int &in_N,
    const double &dt3d )
: N( in_N ), h( dt3d/static_cast<double>(N) ),
  absTol( 1.0e-8 ), relTol( 1.0e-5 )
{
  // Now read the lpn files for num_ebc, Rd, C, and Rp
  SYS_T::file_check( lpn_filename ); // make sure the file is on the disk

  std::ifstream reader;
  reader.open( lpn_filename.c_str(), std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;
  std::string bc_type;

  // The first non-commented line should be
  // RCR num_ebc
  while( std::getline(reader, sline) )
  {
    if( sline[0] != '#' && !sline.empty() )
    {
      sstrm.str(sline);
      sstrm >> bc_type;
      sstrm >> num_ebc;
      sstrm.clear();
      break;
    }
  }

  // Check the file's bc_type matches RCR
  if( bc_type.compare("RCR") ==0 
      || bc_type.compare("rcr") == 0 
      || bc_type.compare("Rcr") == 0 )
  {
    Rd.resize( num_ebc ); C.resize( num_ebc ); 
    Rp.resize( num_ebc ); Pd.resize( num_ebc );
    Q0.resize( num_ebc ); Pi0.resize( num_ebc );
  }
  else SYS_T::print_fatal("Error: the outflow model in %s does not match GenBC_Resistance.\n", lpn_filename.c_str());

  // Read files for each ebc to set the values of Rp, C, and Rd
  int counter = 0;
  while( std::getline(reader, sline) )
  {
    if( sline[0] != '#' && !sline.empty() )
    {
      sstrm.str( sline );
      int face_id;
      sstrm >> face_id;

      // Make sure the face_id, the first column in the file are listed
      // from 0 to ebc_id - 1
      if(face_id != counter) SYS_T::print_fatal("Error: GenBC_RCR the input file %s has wrong format in the face id column (the first column). \n", lpn_filename.c_str());

      sstrm >> Rp[ counter ];
      sstrm >> C[ counter ];
      sstrm >> Rd[ counter ];
      sstrm >> Pd[ counter ];

      sstrm.clear();
      counter += 1;
    }
  }

  if(counter != num_ebc ) SYS_T::print_fatal("Error: GenBC_RCR the input file %s does not contain complete data for outlet faces. \n", lpn_filename.c_str());

  reader.close();

  SYS_T::commPrint( "===> GenBC_RCR data are read in from %s.\n", lpn_filename.c_str() );

  // Set a zero initial value. They should be reset based on the initial
  // 3D solutions.
  for(int ii=0; ii<num_ebc; ++ii)
  {
    Q0[ii] = 0.0; Pi0[ii] = 0.0;

    // Make sure C and Rd are nonzero  
    SYS_T::print_fatal_if(C[ii]==0.0, "Error: GenBC_RCR C cannot be zero.\n");
    SYS_T::print_fatal_if(Rd[ii]==0.0, "Error: GenBC_RCR Rd cannot be zero.\n");
  }
}

void GenBC_RCR::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint( "===> RCR model: N = %d, h = %e, num_ebc = %d \n", N, h, num_ebc );

  for(int ii=0; ii<num_ebc; ++ii)
    SYS_T::commPrint( "     ebcid = %d, Rp = %e, C = %e, Rd =%e, Pd = %e \n", ii, Rp[ii], C[ii], Rd[ii], Pd[ii] );
  
  SYS_T::print_sep_line();
}


double GenBC_RCR::get_m( const int &ii, const double &in_dot_Q,
   const double &in_Q ) const
{
  double diff = std::abs(in_Q) * relTol;

  if( diff < absTol ) diff = absTol;  

  const double left  = get_P(ii, in_dot_Q, in_Q + 0.5 * diff);
  const double right = get_P(ii, in_dot_Q, in_Q - 0.5 * diff);

  return (left - right) / diff;
}


double GenBC_RCR::get_P( const int &ii, const double &in_dot_Q,
   const double &in_Q, const double &time ) const
{
  const double fac13 = 1.0 / 3.0;
  const double fac23 = 2.0 / 3.0;
  const double fac18 = 1.0 / 8.0;
  const double fac38 = 3.0 / 8.0;

  double pi_m = Pi0[ii]; // Pi_m
  
  // in_Q gives Q_N = Q_n+1, and Q0[ii] gives Q_0 = Q_n
  for(int mm=0; mm<N; ++mm)
  {
    const double Q_m = Q0[ii] + static_cast<double>(mm) * ( in_Q - Q0[ii] ) / static_cast<double>(N);
    
    const double Q_mp1 = Q0[ii] + static_cast<double>(mm+1) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

    const double K1 = F(ii, pi_m, Q_m );

    const double K2 = F(ii, pi_m + fac13 * K1 * h, fac23 * Q_m + fac13 * Q_mp1 );

    const double K3 = F(ii, pi_m - fac13 * K1 * h + K2 * h, fac13 * Q_m + fac23 * Q_mp1);
  
    const double K4 = F(ii, pi_m + K1 * h - K2 * h + K3 * h, Q_mp1);

    pi_m = pi_m + fac18 * K1 * h + fac38 * K2 * h + fac38 * K3 * h + fac18 * K4 * h; 
  }

  return pi_m + Rp[ii] * in_Q + Pd[ii];
}


void GenBC_RCR::reset_initial_sol( const int &ii, const double &in_Q_0,
    const double &in_P_0 )
{
  Q0[ii]  = in_Q_0;

  Pi0[ii] = in_P_0 - in_Q_0 * Rp[ii] - Pd[ii];
}

// EOF
