#include "GenBC_Resistance.hpp"

GenBC_Resistance::GenBC_Resistance( const char * const &lpn_filename )
{
  // Now read the values of resis_R and resis_Pd from disk file lpn_filename
  std::string temp_name( lpn_filename );
  SYS_T::file_check( temp_name ); // make sure the file is on the disk

  std::ifstream reader;
  reader.open( lpn_filename, std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;

  std::string bc_type;

  // The first non-commented line should contain
  // bc_type [Resistance, RCR] number-of-ebc
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
  if( bc_type.compare("Resistance") == 0 
      || bc_type.compare("resistance") == 0 
      || bc_type.compare("RESISTANCE") == 0 )
  {
    resis.resize( num_ebc ); pres_offset.resize( num_ebc );
    Q0.resize( num_ebc ); P0.resize( num_ebc );
  }
  else SYS_T::print_fatal("Error: the outflow model in %s does not match GenBC_Resistance.\n", lpn_filename);
 
  // Read files for each ebc to define the parameters for LPN
  int counter = 0;
  while( std::getline(reader, sline) )
  {
    if( sline[0] != '#' && !sline.empty() )
    {
      sstrm.str( sline );
      int face_id;
      sstrm >> face_id;
      
      if(face_id != counter) SYS_T::print_fatal("Error: GenBC_Resistance the input file %s has wrong format in the face id column (the first column). \n", lpn_filename);
      
      sstrm >> resis[ counter ];
      sstrm >> pres_offset[ counter ];

      sstrm.clear();
      
      counter += 1;
    }
  }

  if(counter != num_ebc ) SYS_T::print_fatal("Error: GenBC_Resistance the input file %s does not contain complete data for outlet faces. \n", lpn_filename);

  reader.close();

  SYS_T::commPrint( "===> GenBC_Resistance data are read in from %s.\n", lpn_filename );

  // Set zero initial value.
  for(int ii=0; ii<num_ebc; ++ii)
  {
    Q0[ii] = 0.0;
    P0[ii] = 0.0;
  }
}


GenBC_Resistance::~GenBC_Resistance()
{}


void GenBC_Resistance::print_info() const
{
  SYS_T::commPrint("     GenBC_Resistance : \n");
  for(int ii=0; ii<num_ebc; ++ii)
    SYS_T::commPrint( "     ebcid = %d, R = %e, Pd = %e \n", ii, resis[ii], pres_offset[ii] );
}

// EOF
