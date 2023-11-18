#include "GenBC_Tools.hpp"

int GENBC_T::get_genbc_file_type( const std::string &lpn_filename )
{
  // open the file
  SYS_T::file_check( lpn_filename );

  std::ifstream reader;
  reader.open( lpn_filename.c_str(), std::ifstream::in );

  // Read the fist non-comment line of the file
  std::istringstream sstrm;
  std::string sline, bc_type;

  while( std::getline(reader, sline) )
  {
    if( sline[0] != '#' && !sline.empty() )
    {
      sstrm.str(sline);
      sstrm >> bc_type;
      sstrm.clear();
      break;
    }
  }

  reader.close();

  // Now compare and return the corresponding type value
  if( bc_type.compare("Resistance") == 0
      || bc_type.compare("resistance") == 0
      || bc_type.compare("RESISTANCE") == 0 )
  {
    return 1;
  }
  else if( bc_type.compare("RCR") ==0
      || bc_type.compare("rcr") == 0
      || bc_type.compare("Rcr") == 0 )
  {
    return 2;
  }
  else if( bc_type.compare("Inductance") ==0
      || bc_type.compare("inductance") == 0
      || bc_type.compare("INDUCTANCE") == 0 )
  {
    return 3;
  }
  else if( bc_type.compare("Coronary") ==0
      || bc_type.compare("coronary") == 0
      || bc_type.compare("CORONARY") == 0 )
  {
    return 4;
  }
  else if( bc_type.compare("Pressure") ==0
      || bc_type.compare("pressure") == 0
      || bc_type.compare("PRESSURE") == 0 )
  {
    return 5;
  }
  else if( bc_type.compare("Absorbing") ==0
      || bc_type.compare("absorbing") == 0
      || bc_type.compare("ABSORBING") == 0 )
  {
    return 6;
  }
  else
  {
    return 0;
  }
}

void GENBC_T::set_pchip( const int &np, const std::vector<double> &xp, 
    const std::vector<double> &fp, std::vector<double> &dp )
{
  SYS_T::print_fatal_if( np < 2, "Error: GenBC_Tools set_pchip: Number of evaluation points is less than 2. \n");

  for (int ii=1; ii<np; ++ii)
    SYS_T::print_fatal_if(xp[ii] <= xp[ii-1], "Error: GenBC_Tools set_pchip: xp array not strictly increasing. \n");

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
  }
  else
  {
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
    if( sign_test ( dp[0], del1 ) <= 0.0 ) dp[0] = 0.0;
    else if( sign_test ( del1, del2 ) < 0.0 )
    {
      // Need do this check only if monotonicity switches.
      dmax = 3.0 * del1;

      if ( fabs ( dmax ) < fabs ( dp[0] ) ) dp[0] = dmax;
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

      temp = sign_test ( del1, del2 );

      if( temp < 0.0 )
      {
        ierr = ierr + 1;
        dsave = del2;
      }
      else if( temp == 0.0 )
      {
        // Count number of changes in direction of monotonicity.
        if ( del2 != 0.0 )
        {
          if ( sign_test ( dsave, del2 ) < 0.0 ) ierr = ierr + 1;
          
          dsave = del2;
        }
      }
      else
      {
        // Use Brodlie modification of Butland formula.
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

    if( sign_test ( dp[np-1], del2 ) <= 0.0 )
      dp[np-1] = 0.0;
    else if( sign_test ( del1, del2 ) < 0.0 )
    {
      //
      //  Need do this check only if monotonicity switches.
      //
      dmax = 3.0 * del2;

      if( fabs ( dmax ) < fabs ( dp[np-1] ) )
        dp[np-1] = dmax;
    }
  } // End of If-else statement for np
}

double GENBC_T::sign_test( const double &arg1, const double &arg2 ) 
{
  double value = 0.0;

  if ( arg1 == 0.0 ) value = 0.0;
  else if ( arg1 < 0.0 )
  {
    if ( arg2 < 0.0 ) value = 1.0;
    else if ( arg2 == 0.0 ) value = 0.0;
    else value = -1.0;
  }
  else
  {
    if ( arg2 < 0.0 ) value = -1.0;
    else if ( arg2 == 0.0 ) value = 0.0;
    else value = 1.0;
  }

  return value;
}

void GENBC_T::get_cubic_hermite( const double &x1, const double &x2, 
    const double &f1, const double &f2, const double &d1, const double &d2, 
    const int &ne, const std::vector<double> &xe, std::vector<double> &fe ) 
{
  SYS_T::print_fatal_if(ne<1, "Error: GenBC_Tools get_cubic_hermite: Number of evaluation points is less than 1 \n");

  const double hh = x2 - x1;

  SYS_T::print_fatal_if(hh==0.0, "Error: GenBC_Tools get_cubic_hermite: The interval [x1,x2] is of zero length \n");
  
  const double c1 = hh * d1;

  const double c2 = -3.0 * f1 - 2.0 * hh * d1 + 3.0 * f2 - hh * d2;

  const double c3 = 2.0 * f1 + hh * d1 - 2.0 * f2 + hh * d2;

  // Evaluation loop.
  for (int ii=0; ii<ne; ++ii)
  {
    const double tt = (xe[ii] - x1) / hh;
    fe[ii] = f1 + tt * ( c1 + tt * ( c2 + tt * c3 ) ) ;
  }
}

void GENBC_T::get_cubic_hermite_der( const double &x1, const double &x2, 
    const double &f1, const double &f2, const double &d1, const double &d2, 
    const int &ne, const std::vector<double> &xe, std::vector<double> &de ) 
{
  SYS_T::print_fatal_if(ne<1, "Error: GenBC_Tools get_cubic_hermite_der: Number of evaluation points is less than 1. \n");

  const double hh = x2 - x1;

  SYS_T::print_fatal_if(hh==0.0, "Error: GenBC_Tools get_cubic_hermite_der: The interval [x1,x2] is of zero length. \n");

  const double c2 =-6.0 * f1 / hh - 4.0 * d1 + 6.0 * f2 / hh - 2.0 * d2;

  const double c3 = 6.0 * f1 / hh + 3.0 * d1 - 6.0 * f2 / hh + 3.0 * d2;

  // Evaluation loop.
  for (int ii=0; ii<ne; ++ii)
  {
    const double tt = (xe[ii] - x1)/hh;
    de[ii] = tt * ( c2 + tt * c3 ) + d1;
  }
}

// EOF
