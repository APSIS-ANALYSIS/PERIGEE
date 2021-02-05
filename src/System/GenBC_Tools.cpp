#include "GenBC_Tools.hpp"

int GENBC_T::get_genbc_file_type( const char * const &lpn_filename )
{
  // open the file
  const std::string temp_name( lpn_filename );
  SYS_T::file_check( temp_name );

  std::ifstream reader;
  reader.open( lpn_filename, std::ifstream::in );

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
  else
  {
    return 0;
  }
}

// EOF
