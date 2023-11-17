#include "GenBC_Inductance.hpp"

GenBC_Inductance::GenBC_Inductance( const std::string &lpn_filename )
{
  // Now read the values of induct_L and induct_p from disk file lpn_filename
  SYS_T::file_check( lpn_filename ); // make sure the file is on the disk

  std::ifstream reader;
  reader.open( lpn_filename.c_str(), std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;

  std::string bc_type;

  // The first non-commented line should contain
  // bc_type [Resistance, RCR, Inductance] number-of-ebc
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

  // Allocate the distal pressure and resistance vector
  if( bc_type.compare("Inductance") == 0 
      || bc_type.compare("inductance") == 0 
      || bc_type.compare("INDUCTANCE") == 0 )
  {
    induct.resize( num_ebc ); pres_offset.resize( num_ebc );
    Q0.resize( num_ebc ); P0.resize( num_ebc );
  }
  else SYS_T::print_fatal("Error: the outflow model in %s does not match GenBC_Inductance.\n", lpn_filename.c_str());
 
  // Read files for each ebc to define the parameters for LPN
  int counter = 0;
  while( std::getline(reader, sline) )
  {
    if( sline[0] != '#' && !sline.empty() )
    {
      sstrm.str( sline );
      int face_id;
      sstrm >> face_id;
      
      if(face_id != counter) SYS_T::print_fatal("Error: GenBC_Inductance the input file %s has wrong format in the face id column (the first column). \n", lpn_filename.c_str());
      
      sstrm >> induct[ counter ];
      sstrm >> pres_offset[ counter ];

      sstrm.clear();
      
      counter += 1;
    }
  }

  if(counter != num_ebc ) SYS_T::print_fatal("Error: GenBC_Inductance the input file %s does not contain complete data for outlet faces. \n", lpn_filename.c_str());

  reader.close();

  SYS_T::commPrint( "===> GenBC_Inductance data are read in from %s.\n", lpn_filename.c_str() );

  // Set zero initial value.
  for(int ii=0; ii<num_ebc; ++ii)
  {
    Q0[ii] = 0.0;
    P0[ii] = 0.0;
  }
}

void GenBC_Inductance::print_info() const
{
  SYS_T::commPrint("===> GenBC_Inductance : \n");
  for(int ii=0; ii<num_ebc; ++ii)
    SYS_T::commPrint( "     ebcid = %d, L = %e, p = %e \n", ii, induct[ii], pres_offset[ii] );
}

// EOF
