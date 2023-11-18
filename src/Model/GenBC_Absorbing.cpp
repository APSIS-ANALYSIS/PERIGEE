#include "GenBC_Absorbing.hpp"

GenBC_Absorbing::GenBC_Absorbing( const std::string &lpn_filename,
    const double &solid_E, const double &solid_nu )
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
    para_beta.resize(num_ebc); initial_outlet_area.resize(num_ebc);
  }
  else
    SYS_T::print_fatal( "GenBC_Absorbing Error: BC type in %s should be Absorbing.\n", lpn_filename.c_str() );

  // Read in initial radius and initial thickness per ebc
  for(int ebc_id=0; ebc_id<num_ebc; ++ebc_id)
  { 
    double radius {0}, thickness {0};
    while( std::getline(reader, sline) )
    {
      // face_id num_of_mode w period
      if( sline[0] !='#' && !sline.empty() )
      {
        sstrm.str(sline);
        int face_id;
        sstrm >> face_id;
        SYS_T::print_fatal_if( face_id != ebc_id, "GenBC_Pressure Error: ebc in %s should be listed in ascending order.\n", lpn_filename.c_str() );

        sstrm >> radius;
        SYS_T::print_fatal_if(radius <=0, "GenBC_Absorbing Error: Radius of ebc_id %d should be greater than 0.\n", ebc_id);

        sstrm >> thickness;
        SYS_T::print_fatal_if(radius <=0, "GenBC_Absorbing Error: Thickness of ebc_id %d should be greater than 0.\n", ebc_id);

        sstrm.clear();
        break;
      }
    }

    // calculate beta and outlet area
    para_beta[ebc_id] = thickness * (solid_E / (1 - solid_nu * solid_nu)) / (radius * radius);

    initial_outlet_area[ebc_id] = MATH_T::PI * radius * radius;

    current_outlet_area[ebc_id] = initial_outlet_area[ebc_id];
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
    SYS_T::commPrint(" beta = %e, Initial surface area =%e \n", para_beta[ii], initial_outlet_area[ii]);
  }
}

double GenBC_Absorbing::get_P( const int &ii, const double &dot_Q, const double &Q,
    const double &time) const
{
  return para_beta[ii] * (std::sqrt(current_outlet_area[ii]) - std::sqrt(initial_outlet_area[ii])) / MATH_T::PI;
}

// EOF