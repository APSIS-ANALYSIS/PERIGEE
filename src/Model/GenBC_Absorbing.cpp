#include "GenBC_Absorbing.hpp"

GenBC_Absorbing::GenBC_Absorbing( const std::string &lpn_filename,
    const double &solid_E, const double &solid_nu, const double &in_fl_density )
{
  SYS_T::file_check( lpn_filename ); // make sure the file is on the disk

  std::ifstream reader;
  reader.open( lpn_filename.c_str(), std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;
  std::string bc_type;

  // The first non-commented, non-empty line should be
  // Absorbing num_ebc
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

  if( bc_type.compare("Absorbing") == 0 || bc_type.compare("ABSORBING") == 0 || bc_type.compare("absorbing") == 0)
  {
    para_1.resize(num_ebc); para_2.resize(num_ebc);
    P_ref.resize(num_ebc); P0.resize(num_ebc); P.resize(num_ebc);
  }
  else
    SYS_T::print_fatal( "GenBC_Absorbing Error: BC type in %s should be Absorbing.\n", lpn_filename.c_str() );

  // Read in initial radius and initial thickness per ebc
  for(int ebc_id=0; ebc_id<num_ebc; ++ebc_id)
  { 
    double A0 {0}, A_bar {0}, thickness {0};
    while( std::getline(reader, sline) )
    {
      // face_id num_of_mode w period
      if( sline[0] !='#' && !sline.empty() )
      {
        sstrm.str(sline);
        int face_id;
        sstrm >> face_id;
        SYS_T::print_fatal_if( face_id != ebc_id, "GenBC_Pressure Error: ebc in %s should be listed in ascending order.\n", lpn_filename.c_str() );

        sstrm >> thickness;
        SYS_T::print_fatal_if(thickness <= 0, "GenBC_Absorbing Error: Thickness of ebc_id %d should be greater than 0.\n", ebc_id);

        sstrm >> A0;
        SYS_T::print_fatal_if(A0 <= 0, "GenBC_Absorbing Error: Initial surface area of ebc_id %d should be greater than 0.\n", ebc_id);

        sstrm >> A_bar;
        SYS_T::print_fatal_if(A_bar <= 0, "GenBC_Absorbing Error: Steady-state surface area of ebc_id %d should be greater than 0.\n", ebc_id);

        sstrm >> P_ref[ebc_id];

        sstrm >> P0[ebc_id];

        sstrm.clear();
        break;
      }
    }

    // calculate parameters
    para_1[ebc_id] = thickness * (solid_E / (1 - solid_nu * solid_nu)) * MATH_T::PI / std::sqrt(A0);

    para_2[ebc_id] = std::sqrt(in_fl_density / 8) * (1.0 / A_bar);
  }

  // Finish reading the file and close it
  reader.close();
}

void GenBC_Absorbing::print_info() const
{
  SYS_T::commPrint("===> GenBC_Absorbing: \n");

  for(int ii=0; ii<num_ebc; ++ii)
  {
    SYS_T::commPrint("     ebc_id = %d", ii);
    SYS_T::commPrint(" Parameter 1 = %e, Parameter 2 =  %e, P_ref = %e\n",
      para_1[ii], para_2[ii], P_ref[ii]);
  }
}

double GenBC_Absorbing::set_P( const int &ii, const double &Q0) const
{
  const double temp = (para_2[ii] * Q0 + std::sqrt(para_1[ii]));
  return (temp * temp - para_1[ii] + P_ref[ii]);
}

// EOF