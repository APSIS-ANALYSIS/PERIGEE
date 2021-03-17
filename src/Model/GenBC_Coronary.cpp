#include "GenBC_Coronary.hpp"

GenBC_Coronary::GenBC_Coronary( const char * const &lpn_filename, 
    const int &in_N, const double &dt3d )
: num_odes(2), N( in_N ), h( dt3d/static_cast<double>(N) ), 
  absTol( 1.0e-8 ), relTol( 1.0e-5 )
{
  // Now read the lpn input file for num_ebc and coronary model 
  // parameters (Ra, Ca, Ra_micro, Cim, Rv, Pd and Pim)
  std::string temp_name( lpn_filename );
  SYS_T::file_check( temp_name ); // make sure the file is on the disk

  std::ifstream reader;
  reader.open( lpn_filename, std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;
  std::string bc_type;

  // The first non-commented line should be Coronary num_ebc
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

  // Check the file's bc_type matches Coronary
  if( bc_type.compare("Coronary") ==0
      || bc_type.compare("CORONARY") == 0
      || bc_type.compare("coronary") == 0 )
  {
    // Allocate the containers for model parameters
    Ra.resize( num_ebc ); 
    Ca.resize( num_ebc );
    Ra_micro.resize( num_ebc ); 
    Cim.resize( num_ebc );
    Rv.resize( num_ebc ); 
    Pd.resize( num_ebc );
    alpha_Pim.resize( num_ebc );
    Q0.resize( num_ebc );
    Pi0.resize( num_ebc ); 
    num_Pim_data.resize( num_ebc );
    Time_data.resize( num_ebc );
    Pim_data.resize( num_ebc );
    der_Pim_data.resize( num_ebc );
    prev_0D_sol.resize( num_ebc );
    dPimdt_k1.resize( num_ebc );
    dPimdt_k2.resize( num_ebc );
    dPimdt_k3.resize( num_ebc );

    for( int ii =0; ii<num_ebc; ++ii )
    {
      prev_0D_sol[ii].resize( num_odes );
      Pi0[ii].resize( num_odes );

      dPimdt_k1[ii].resize( N+1 );
      dPimdt_k2[ii].resize( N );
      dPimdt_k3[ii].resize( N );
    }
  }
  else SYS_T::print_fatal("Error: the outflow model in %s does not match GenBC_Coronary.\n", lpn_filename);

  // Read files for each ebc to set the values of Ra, Ca, Ra_micro, 
  // Cim, Rv, Pd, num_Pim_data, and alpha_Pim
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
      SYS_T::print_fatal_if( face_id != counter, 
          "Error: GenBC_Coronary the input file %s has wrong format in the face id column (the first column). \n", lpn_filename );

      sstrm >> Ra[ counter ];
      sstrm >> Ca[ counter ];
      sstrm >> Ra_micro[ counter ];
      sstrm >> Cim[ counter ];
      sstrm >> Rv[ counter ];
      sstrm >> Pd[ counter ];
      sstrm >> num_Pim_data[ counter ];
      sstrm >> alpha_Pim[ counter ];

      SYS_T::print_fatal_if( num_Pim_data[counter] <= 2 && num_Pim_data[counter] != 0, 
          "Error: number of Pim data needs to be 0 for RCR or greater than 2 for coronary BC. \n" );

      // Resize the IntraMyocardial data 
      if( num_Pim_data[counter] > 0 ) 
      {
        is_RCR = false;

        const int data_size = num_Pim_data[counter];

        Time_data[counter].resize( data_size );
        Pim_data[counter].resize( data_size );
        der_Pim_data[counter].resize( data_size );

        sstrm.clear();

        for(int ii =0; ii < num_Pim_data[counter]; ++ii)
        {
          getline(reader, sline);
          sstrm.str( sline );
          sstrm>>Time_data[counter][ii];
          sstrm>>Pim_data[counter][ii];

          Pim_data[counter][ii] = Pim_data[counter][ii] * alpha_Pim[counter];
          sstrm.clear();
        }

        GENBC_T::set_pchip( num_Pim_data[counter], Time_data[counter], Pim_data[counter], der_Pim_data[counter] );
        get_dPim_dt( counter, 0.0, N * h );

        SYS_T::print_fatal_if( Time_data[counter][0]>0.0, "Error: Pim data does not start from 0.\n" );
      }
      else
      {
        is_RCR = true;

        Time_data[counter].clear();
        Pim_data[counter].clear();
        der_Pim_data[counter].clear();
      }
      
      sstrm.clear();
      counter += 1;
    }
  }

  SYS_T::print_fatal_if( counter != num_ebc, "Error: GenBC_Coronary the input file %s does not contain complete data for outlet faces.\n", lpn_filename );

  reader.close();

  SYS_T::commPrint( "===> GenBC_Coronary data are read in from %s.\n", lpn_filename );

  // Set zero initial values. 
  // They will be reset based on the 3D solutions.
  for(int ii=0; ii<num_ebc; ++ii)
  {
    Q0[ii] = 0.0;
    for(int jj=0; jj<num_odes; ++jj)
    {
      Pi0[ii][jj] = 0.0;
      prev_0D_sol[ii][jj] = 0.0;
    }  

    // Make sure C and R are nonzero
    SYS_T::print_fatal_if( Ca[ii] == 0.0, "Error: GenBC_Coronary Ca cannot be zero.\n" );
    SYS_T::print_fatal_if( Cim[ii] == 0.0, "Error: GenBC_Coronary Cim cannot be zero.\n" );
    SYS_T::print_fatal_if( Ra_micro[ii] == 0.0, "Error: GenBC_Coronary Ra_micro cannot be zero.\n" );
    SYS_T::print_fatal_if( Rv[ii] == 0.0, "Error: GenBC_Coronary Rv cannot be zero.\n" );
  }
}

GenBC_Coronary::~GenBC_Coronary()
{}

void GenBC_Coronary::print_info() const
{
  SYS_T::commPrint( "     Coronary model: N = %d, h = %e, num_ebc = %d \n", N, h, num_ebc );

  for(int ii=0; ii<num_ebc; ++ii)
    SYS_T::commPrint( "     ebcid = %d, Ra = %e, Ca = %e, Ra_micro = %e,Rim = %e, Rv = %e, Pd = %e \n", 
        ii, Ra[ii], Ca[ii],Ra_micro[ii],Cim[ii], Rv[ii], Pd[ii] );
}

double GenBC_Coronary::get_m( const int &ii, const double &in_dot_Q, const double &in_Q ) const
{
  double diff = std::abs(in_Q) * relTol;

  if( diff < absTol ) diff = absTol;

  const double left  = get_P(ii, in_dot_Q, in_Q + 0.5 * diff);
  const double right = get_P(ii, in_dot_Q, in_Q - 0.5 * diff);

  return (left - right) / diff;
}

double GenBC_Coronary::get_P( const int &ii, const double &in_dot_Q,
    const double &in_Q ) const
{
  double output_P = 0.0;

  const double fac13 = 1.0 / 3.0;
  const double fac23 = 2.0 / 3.0;
  const double fac18 = 1.0 / 8.0;
  const double fac38 = 3.0 / 8.0;

  if( is_RCR )
  {
    // Each RCR face is governed by an ODE only.
    double pi_m = Pi0[ii][0]; // Pi_m

    // in_Q gives Q_N = Q_n+1, and Q0[ii] gives Q_0 = Q_n. 
    // do Runge-Kutta 4 with the 3/8 rule
    for(int mm=0; mm<N; ++mm)
    {
      const double Q_m = Q0[ii] + static_cast<double>(mm) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

      const double Q_mp1 = Q0[ii] + static_cast<double>(mm+1) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

      const double K1 = F_RCR(ii, pi_m, Q_m);

      const double K2 = F_RCR(ii, pi_m + fac13 * K1 * h, fac23 * Q_m + fac13 * Q_mp1);

      const double K3 = F_RCR(ii, pi_m - fac13 * K1 * h + K2 * h, fac13 * Q_m + fac23 * Q_mp1);

      const double K4 = F_RCR(ii, pi_m + K1 * h - K2 * h + K3 * h, Q_mp1);

      pi_m = pi_m + fac18 * K1 * h + fac38 * K2 * h + fac38 * K3 * h + fac18 * K4 * h;
    }

    prev_0D_sol[ii][0] = pi_m;

    output_P =  pi_m + Ra[ii] * in_Q ; 
  }
  else
  {
    // Here, we know it is a coronary face.
    // Each coronary face is governed by num_odes = 2 ODEs.  
    // initial pressures at Ca and Cim
    double pi_m[num_odes];

    // auxiliary variables for RK4 
    double K1[num_odes], K2[num_odes], K3[num_odes], K4[num_odes];
    double pi_tmp[num_odes];

    for(int jj=0; jj<num_odes; ++jj)
    {
      pi_m[jj] = Pi0[ii][jj];
      K1[jj] = 0.0;
      K2[jj] = 0.0;
      K3[jj] = 0.0;
      K4[jj] = 0.0;
      pi_tmp[jj] = 0.0;
    }

    // in_Q gives Q_N = Q_n+1, and Q0[ii] gives Q_0 = Q_n
    // do Runge-Kutta 4 with the 3/8 rule
    for(int mm=0; mm<N; ++mm)
    {
      const double Q_m = Q0[ii] + static_cast<double>(mm) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

      const double Q_mp1 = Q0[ii] + static_cast<double>(mm+1) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

      F_coronary(ii, pi_m, Q_m, dPimdt_k1[ii][mm], K1);

      for(int jj=0; jj<num_odes; ++jj)
        pi_tmp[jj] = pi_m[jj] + fac13 * K1[jj] * h;

      F_coronary(ii, pi_tmp , fac23 * Q_m + fac13 * Q_mp1, dPimdt_k2[ii][mm], K2);

      for(int jj=0; jj<num_odes; ++jj)
        pi_tmp[jj] = pi_m[jj] - fac13*K1[jj] * h + K2[jj] * h;

      F_coronary(ii, pi_tmp , fac13 * Q_m + fac23 * Q_mp1, dPimdt_k3[ii][mm], K3);

      for(int jj=0; jj<num_odes; ++jj)
        pi_tmp[jj] = pi_m[jj] + K1[jj] * h - K2[jj] * h + K3[jj] * h;

      F_coronary(ii, pi_tmp, Q_mp1, dPimdt_k1[ii][mm+1], K4);

      for(int jj=0; jj<2; ++jj)
        pi_m[jj] = pi_m[jj] + fac18 * K1[jj] * h + fac38 * K2[jj] * h + fac38 * K3[jj] * h + fac18 * K4[jj] * h;
    }

    // Make a copy of the ODE solutions. 
    // prev_0D_sol will reset initial values Pi0 for future time integration (t=n+1-> t=n+2) 
    for(int jj=0; jj<num_odes; ++jj) prev_0D_sol[ii][jj]=pi_m[jj];

    output_P = pi_m[0] + Ra[ii] * in_Q;
  }
  return output_P;
}

double GenBC_Coronary::get_P0( const int &ii ) const
{
  return Q0[ii] * Ra[ii] + Pi0[ii][0];
}

void GenBC_Coronary::reset_initial_sol( const int &ii, const double &in_Q_0,
    const double &in_P_0, const double &curr_time )
{
  Q0[ii] = in_Q_0;

  // Use the last 0D solition as initial solutions for the next time integration.
  for(int jj=0; jj<num_odes; ++jj) Pi0[ii][jj] = prev_0D_sol[ii][jj];

  // Precalculate dPimdt values needed for integrating Coronary ODEs.
  if( num_Pim_data[ii]>0 ) get_dPim_dt(ii, curr_time, curr_time + N * h);
}

void GenBC_Coronary::F_coronary( const int &ii, const double * const &pi, const double &q, 
    const double &dPimdt, double * const &K ) const
{
  // The Coronary LPM consists of two ODEs. 
  K[0]=(q-(pi[0]-pi[1])/Ra_micro[ii])/Ca[ii];
  K[1]=((pi[0]-pi[1])/Ra_micro[ii]-(pi[1]-Pd[ii])/(Rv[ii]))/Cim[ii]+dPimdt;
}

double GenBC_Coronary::F_RCR( const int &ii, const double &pi, const double &q ) const
{
  // This is an RCR ODE, Ra_micro and Ca are used to store distal resistance and capacitance
  return (q-(pi-Pd[ii])/Ra_micro[ii])/Ca[ii];
}

void GenBC_Coronary::get_dPim_dt( const int &ii, const double &time_start, const double &time_end )
{
  double tend_mod = fmod( time_end, Time_data[ii][num_Pim_data[ii]-1] );
  if ( tend_mod<absTol ) tend_mod =  Time_data[ii][num_Pim_data[ii]-1];

  // Find the interval in Pim that covers current integration time. 
  double x1,x2,f1,f2,d1,d2;
  for(int mm=1; mm<num_Pim_data[ii];++mm)
  {
    if ( tend_mod <= Time_data[ii][mm] )
    {
      x1 = Time_data[ii][mm-1];
      f1 = Pim_data[ii][mm-1];
      d1 = der_Pim_data[ii][mm-1];
      x2 = Time_data[ii][mm];
      f2 = Pim_data[ii][mm];
      d2 = der_Pim_data[ii][mm];
      break;
    }
  }

  std::vector<double> xe_1,xe_2,xe_3;
  xe_1.resize( N+1 );
  xe_2.resize( N );
  xe_3.resize( N );

  // Precalculate dPimdt values at time points evaluated by RK4 
  // from time_start to time_end. 
  double tmp = time_start;

  const double fac13 = 1.0 / 3.0;
  const double fac23 = 2.0 / 3.0;

  for(int mm=0;mm<N;++mm)
  {
    xe_1[mm] = fmod(tmp,Time_data[ii][num_Pim_data[ii]-1]);
    xe_2[mm] = fmod(tmp+fac13*h,Time_data[ii][num_Pim_data[ii]-1]);
    xe_3[mm] = fmod(tmp+fac23*h,Time_data[ii][num_Pim_data[ii]-1]);
    tmp = tmp + h;
  }

  xe_1[N] = tend_mod;

  GENBC_T::get_cubic_hermite_der( x1, x2, f1,f2, d1, d2, N+1, xe_1, dPimdt_k1[ii] );

  GENBC_T::get_cubic_hermite_der( x1, x2, f1,f2, d1, d2, N,   xe_2, dPimdt_k2[ii] );

  GENBC_T::get_cubic_hermite_der( x1, x2, f1,f2, d1, d2, N,   xe_3, dPimdt_k3[ii] );
}

// EOF
