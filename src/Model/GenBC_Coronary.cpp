#include "GenBC_Coronary.hpp"

GenBC_Coronary::GenBC_Coronary( const char * const &lpn_filename, 
    const int &in_N, const double &dt3d )
: num_odes(2), N( in_N ), h( dt3d/static_cast<double>(N) ), 
  absTol( 1.0e-8 ), relTol( 1.0e-5 ),
  tstart( 0.0 ), tend( dt3d )
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

      // !! WHY N+1 HERE?  
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

        spline_pchip_set( num_Pim_data[counter], Time_data[counter], Pim_data[counter], der_Pim_data[counter] );
        get_dPim_dt( counter );

        SYS_T::print_fatal_if(Time_data[counter][0]>0.0, "Error: Pim data does not start from 0.\n");
      }
      else
      {
        Time_data[counter].clear();
        Pim_data[counter].clear();
        der_Pim_data[counter].clear();
      }
      
      counter += 1;
    }
  }

  SYS_T::print_fatal_if( counter != num_ebc, "Error: GenBC_Coronary the input file %s does not contain complete data for outlet faces.\n", lpn_filename);

  reader.close();

  SYS_T::commPrint( "===> GenBC_Coronary data are read in from %s.\n", lpn_filename );

  // Set a zero initial values. 
  // They will be reset based on the 3D solutions.
  for(int ii=0; ii<num_ebc; ++ii)
  {
    Q0[ii] = 0.0;
    Pi0[ii][0] = 0.0;
    Pi0[ii][1]=0.0;
    prev_0D_sol[ii][0]=0.0;
    prev_0D_sol[ii][1]=0.0;

    // Make sure C and R are nonzero
    SYS_T::print_fatal_if(Ca[ii]==0.0, "Error: GenBC_Coronary Ca cannot be zero.\n");
    SYS_T::print_fatal_if(Cim[ii]==0.0, "Error: GenBC_Coronary Cim cannot be zero.\n");
    SYS_T::print_fatal_if(Ra_micro[ii]==0.0, "Error: GenBC_Coronary Ra_micro cannot be zero.\n");
    SYS_T::print_fatal_if(Rv[ii]==0.0, "Error: GenBC_Coronary Rv cannot be zero.\n");
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

double GenBC_Coronary::get_m( const int &ii, const double &in_dot_Q,
    const double &in_Q ) const
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
  const double fac13 = 1.0 / 3.0;
  const double fac23 = 2.0 / 3.0;
  const double fac18 = 1.0 / 8.0;
  const double fac38 = 3.0 / 8.0;

  // If no Pim data is provided, it is an RCR face. 
  if( num_Pim_data[ii] == 0 )
  {
    // Each RCR face is governed by an ODE only.
    double pi_m = Pi0[ii][0]; // Pi_m

    // in_Q gives Q_N = Q_n+1, and Q0[ii] gives Q_0 = Q_n. 
    // do Runge-Kutta 4 with the 3/8 rule
    for(int mm=0; mm<N; ++mm)
    {
      const double Q_m = Q0[ii] + static_cast<double>(mm) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

      const double Q_mp1 = Q0[ii] + static_cast<double>(mm+1) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

      const double K1 = F(ii, pi_m, Q_m);

      const double K2 = F(ii, pi_m + fac13 * K1 * h, fac23 * Q_m + fac13 * Q_mp1);

      const double K3 = F(ii, pi_m - fac13 * K1 * h + K2 * h, fac13 * Q_m + fac23 * Q_mp1);

      const double K4 = F(ii, pi_m + K1 * h - K2 * h + K3 * h, Q_mp1);

      pi_m = pi_m + fac18 * K1 * h + fac38 * K2 * h + fac38 * K3 * h + fac18 * K4 * h;
    }

    prev_0D_sol[ii][0] = pi_m;

    return pi_m + Ra[ii] * in_Q ; 
  }

  // Here, we know it is a coronary face.
  // Each coronary face is governed by 2 ODEs.  
  // initial pressures at Ca and Cim
  double pi_m[num_odes] = {Pi0[ii][0], Pi0[ii][1]};  //Pi_m

  // auxiliary variables for RK4 
  double K1[num_odes] = {0.0, 0.0};
  double K2[num_odes] = {0.0, 0.0};
  double K3[num_odes] = {0.0, 0.0};
  double K4[num_odes] = {0.0, 0.0};
  double pi_tmp[num_odes] ={0.0, 0.0};

  // in_Q gives Q_N = Q_n+1, and Q0[ii] gives Q_0 = Q_n
  // do Runge-Kutta 4 with the 3/8 rule
  for(int mm=0; mm<N; ++mm)
  {
    const double Q_m = Q0[ii] + static_cast<double>(mm) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

    const double Q_mp1 = Q0[ii] + static_cast<double>(mm+1) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

    F(ii, pi_m, Q_m, dPimdt_k1[ii][mm], K1);

    for(int jj=0; jj<num_odes; ++jj)
      pi_tmp[jj] = pi_m[jj] + fac13 * K1[jj] * h;

    F(ii, pi_tmp , fac23 * Q_m + fac13 * Q_mp1, dPimdt_k2[ii][mm], K2);

    for(int jj=0; jj<num_odes; ++jj)
      pi_tmp[jj] = pi_m[jj] - fac13*K1[jj] * h + K2[jj] * h;

    F(ii, pi_tmp , fac13 * Q_m + fac23 * Q_mp1, dPimdt_k3[ii][mm], K3);

    for(int jj=0; jj<num_odes; ++jj)
      pi_tmp[jj] = pi_m[jj] + K1[jj] * h - K2[jj] * h + K3[jj] * h;

    F(ii, pi_tmp, Q_mp1, dPimdt_k1[ii][mm+1], K4);

    for(int jj=0; jj<2; ++jj)
      pi_m[jj] = pi_m[jj] + fac18 * K1[jj] * h + fac38 * K2[jj] * h + fac38 * K3[jj] * h + fac18 * K4[jj] * h;
  }

  // Make a copy of the ODE solutions. 
  // prev_0D_sol will reset initial values Pi0 for future time integration (t=n+1-> t=n+2) 
  for(int jj=0; jj<num_odes; ++jj) prev_0D_sol[ii][jj]=pi_m[jj];

  return pi_m[0] + Ra[ii] * in_Q;
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
  Pi0[ii][0] = prev_0D_sol[ii][0];
  Pi0[ii][1] = prev_0D_sol[ii][1];

  // Update tstart and tend only once. 
  if( ii == 0 )
  {
    tstart=curr_time;
    tend=curr_time+N*h;
  }

  // Precalculate dPimdt values needed for integrating Coronary ODEs.
  if( num_Pim_data[ii]>0 ) get_dPim_dt(ii);
}

void GenBC_Coronary::F( const int &ii, const double * const &pi, const double &q, 
    const double &dPimdt, double * const &K ) const
{
  // The Coronary LPM consists of two ODEs. 
  K[0]=(q-(pi[0]-pi[1])/Ra_micro[ii])/Ca[ii];
  K[1]=((pi[0]-pi[1])/Ra_micro[ii]-(pi[1]-Pd[ii])/(Rv[ii]))/Cim[ii]+dPimdt;
}

double GenBC_Coronary::F(const int &ii, const double &pi, const double &q) const
{
  // This is an RCR ODE, Ra_micro and Ca are used to store distal resistance and capacitance
  return (q-(pi-Pd[ii])/Ra_micro[ii])/Ca[ii];
}


void GenBC_Coronary::spline_pchip_set (const int &np, const std::vector<double> &xp, 
    const std::vector<double> &fp, std::vector<double> &dp)
  //=============================================================================
  // This function sets derivatives for a piecewise cubic Hermite interpolant.
  // This function is modified from John Burkardt's C++ version of the original 
  // Fortran program by Fred Fritsch under the GNU LGPL license.
  //
  // Input: np is the number of points
  //        xp is the corresponding (time) point valuesi, with length np
  //        fp is the corresponding (pressure) point values, with length np
  // Output: dp stores the derivative at xp, with length np
  //
  //
  // Reference:
  //
  // Fred Fritsch, Ralph Carlson,
  //    Monotone Piecewise Cubic Interpolation,
  //    SIAM Journal on Numerical Analysis,
  //    Volume 17, Number 2, April 1980, pages 238-246.
  //
  //    Fred Fritsch, Judy Butland,
  //    A Method for Constructing Local Monotone Piecewise
  //    Cubic Interpolants,
  //    SIAM Journal on Scientific and Statistical Computing,
  //    Volume 5, Number 2, 1984, pages 300-304.
  //
{
  SYS_T::print_fatal_if(np<2, "Error: GenBC_Coronary SPLINE_PCHIP_SET: Number of evaluation points is less than 1 \n");

  for (int ii=1; ii<np; ++ii)
    SYS_T::print_fatal_if(xp[ii] <= xp[ii-1], "Error: GenBC_Coronary SPLINE_PCHIP_SET: X array not strictly increasing. \n");

  int ierr = 0;
  int nless1 = np - 1;
  double h1 = xp[1] - xp[0];
  double del1 = ( fp[1] - fp[0] ) / h1;
  double dsave = del1;

  //  Special case np = 2, use linear interpolation.
  if( np == 2 )
  {
    dp[0] = del1;
    dp[np-1] = del1;
    return;
  }

  // Normal case, np >= 3.
  double h2 = xp[2] - xp[1];
  double del2 = ( fp[2] - fp[1] ) / h2;

  double hsum = h1 + h2;
  double w1 = ( h1 + hsum ) / hsum;
  double w2 = -h1 / hsum;
  double dmax;
  double dmin;
  dp[0] = w1 * del1 + w2 * del2;
  // Set dp[0] via non-centered three point formula, adjusted to be shape preserving.
  if( pch_sign_testing ( dp[0], del1 ) <= 0.0 )
  {
    dp[0] = 0.0;
  }
  // Need do this check only if monotonicity switches.
  else if( pch_sign_testing ( del1, del2 ) < 0.0 )
  {
    dmax = 3.0 * del1;

    if ( fabs ( dmax ) < fabs ( dp[0] ) )
    {
      dp[0] = dmax;
    }
  }

  double temp, hsumt3, drat1, drat2;

  // Loop through interior points.
  for(int ii=2; ii<=nless1; ++ii)
  {
    if( 2 < ii )
    {
      h1 = h2;
      h2 = xp[ii] - xp[ii-1];
      hsum = h1 + h2;
      del1 = del2;
      del2 = ( fp[ii] - fp[ii-1] ) / h2;
    }

    // Set dp[ii-1]=0 unless data are strictly monotonic.
    dp[ii-1] = 0.0;

    temp = pch_sign_testing ( del1, del2 );

    if( temp < 0.0 )
    {
      ierr = ierr + 1;
      dsave = del2;
    }
    // Count number of changes in direction of monotonicity.
    else if( temp == 0.0 )
    {
      if ( del2 != 0.0 )
      {
        if ( pch_sign_testing ( dsave, del2 ) < 0.0 )
        {
          ierr = ierr + 1;
        }
        dsave = del2;
      }
    }
    //  Use Brodlie modification of Butland formula.
    else
    {
      hsumt3 = 3.0 * hsum;
      w1 = ( hsum + h1 ) / hsumt3;
      w2 = ( hsum + h2 ) / hsumt3;
      dmax = fmax ( fabs ( del1 ), fabs ( del2 ) );
      dmin = fmin ( fabs ( del1 ), fabs ( del2 ) );
      drat1 = del1 / dmax;
      drat2 = del2 / dmax;
      dp[ii-1] = dmin / ( w1 * drat1 + w2 * drat2 );
    }
  }

  // Set dp[np-1] via non-centered three point formula, adjusted to be shape preserving.
  w1 = -h2 / hsum;
  w2 = ( h2 + hsum ) / hsum;
  dp[np-1] = w1 * del1 + w2 * del2;

  if( pch_sign_testing ( dp[np-1], del2 ) <= 0.0 )
  {
    dp[np-1] = 0.0;
  }
  else if( pch_sign_testing ( del1, del2 ) < 0.0 )
  {
    //
    //  Need do this check only if monotonicity switches.
    //
    dmax = 3.0 * del2;

    if( fabs ( dmax ) < fabs ( dp[np-1] ) )
    {
      dp[np-1] = dmax;
    }

  }
  return;
}


double GenBC_Coronary::pch_sign_testing ( const double &arg1, const double &arg2 ) const
//=============================================================================
// This function performs a sign test.
// This function is modified from John Burkardt's C++ version of the original 
// Fortran program by Fred Fritsch under the GNU LGPL license.
//
// return -1.0, if arg1 and arg2 are of opposite sign.
// return  0.0, if either argument is zero.
// return +1.0, if arg1 and arg2 are of the same sign.  
//
// The function is to do this without multiplying ARG1 * ARG2, to avoid possible 
// over/underflow problems.
//=============================================================================
{
  double value;

  if ( arg1 == 0.0 )
  {
    value = 0.0;
  }
  else if ( arg1 < 0.0 )
  {
    if ( arg2 < 0.0 )
    {
      value = 1.0;
    }
    else if ( arg2 == 0.0 )
    {
      value = 0.0;
    }
    else if ( 0.0 < arg2 )
    {
      value = -1.0;
    }
  }
  else if ( 0.0 < arg1 )
  {
    if ( arg2 < 0.0 )
    {
      value = -1.0;
    }
    else if ( arg2 == 0.0 )
    {
      value = 0.0;
    }
    else if ( 0.0 < arg2 )
    {
      value = 1.0;
    }
  }

  return value;
}

void GenBC_Coronary::get_dPim_dt(const int &ii)
{
  double x1,x2,f1,f2,d1,d2;

  const double fac13 = 1.0 / 3.0;
  const double fac23 = 2.0 / 3.0;

  double tend_mod = fmod(tend,Time_data[ii][num_Pim_data[ii]-1]);
  // Find the interval in Pim that covers current integration time. 
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
  // from tstart to tend. 
  double tmp=tstart;

  for(int mm=0;mm<N;++mm)
  {
    xe_1[mm] = fmod(tmp,Time_data[ii][num_Pim_data[ii]-1]);
    xe_2[mm] = fmod(tmp+fac13*h,Time_data[ii][num_Pim_data[ii]-1]);
    xe_3[mm] = fmod(tmp+fac23*h,Time_data[ii][num_Pim_data[ii]-1]);
    tmp=tmp+h;
  }

  xe_1[N] = tend_mod;

  cubic_hermite_derivative ( x1, x2, f1,f2, d1, d2, N+1, xe_1, dPimdt_k1[ii] );

  cubic_hermite_derivative ( x1, x2, f1,f2, d1, d2, N, xe_2, dPimdt_k2[ii] );

  cubic_hermite_derivative ( x1, x2, f1,f2, d1, d2, N, xe_3, dPimdt_k3[ii] );
}


void GenBC_Coronary::cubic_hermite_derivative( const double &x1, const double &x2, 
    const double &f1, const double &f2, const double &d1, const double &d2, 
    const int &ne, const std::vector<double> &xe, std::vector<double> &de ) const
//=============================================================================
// A cubic Hermite spline in [x1 x2] with starting and ending points f1 and f2 
// and derivatives d1 and d2 is given by
// f(x)=f1*h00(t)+delta*d1*h10(t)+f2*h01(t)+delta*d2*h11(t),
// where h00(t)=2*t^3-3*t^2+1,h10(t)=t^3-2*t^2+t,h01(t)=-2*t^3+3*t^2,
// h11(t)=t^3-t^2 and delta=x2-x1,t=(x-x1)/delta. 
// This function calculates df(x)/dx for x in [x1 x2].
//=============================================================================    
{

  SYS_T::print_fatal_if(ne<1, "Error: GenBC_Coronary cubic_hermite_derivative: Number of evaluation points is less than 1 \n");

  const double hh = x2 - x1;

  SYS_T::print_fatal_if(hh==0.0, "Error: GenBC_Coronary cubic_hermite_derivative: The interval [X1,X2] is of zero length \n");

  const double c2 =-6.0 * f1 / hh - 4.0 * d1 + 6.0 * f2 / hh - 2.0 * d2;

  const double c3 = 6.0 * f1 / hh + 3.0 * d1 - 6.0 * f2 / hh + 3.0 * d2;

  // Evaluation loop.
  for (int ii=0; ii<ne; ++ii)
  {
    const double tt = (xe[ii] - x1)/hh;
    de[ii] = tt * ( c2 + tt * c3 )+d1 ;
  }
}

// EOF
