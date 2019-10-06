// ==================================================================
// This is a file converter that reads in the rcrt.dat file and 
// convert it to the format for lpn_rcr_input.txt file,
// which lists face-id, Rp, C, Rd, Pd by lines.
//
// Author: Ju Liu
// Date: Sept 13 2019
// ==================================================================
#include "Sys_Tools.hpp"

using std::cout;
using std::cerr;
using std::endl;

int main( int argc, char * argv [] )
{
  std::string sv_rcr_file("rcrt.dat");
  std::string out_rcr_file("lpn_rcr_input.txt");
  int num_outlet = 1;

  PetscMPIInt size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  SYS_T::print_fatal_if(size!=1,"ERROR: converter is a serial routine! \n");

  SYS_T::GetOptionString("-sv_rcr_file", sv_rcr_file);
  SYS_T::GetOptionInt("-num_outlet", num_outlet);

  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -sv_rcr_file: "<<sv_rcr_file<<endl;
  cout<<" -num_outlet: "<<num_outlet<<endl;
  cout<<"names to be written: \n";
  cout<<" -out_rcr_file: "<<out_rcr_file<<endl;
  cout<<"================================"<<endl;

  // Make sure the file exists on disk
  SYS_T::file_check( sv_rcr_file );

  std::ifstream linecounter( sv_rcr_file.c_str(), std::ifstream::in );
  std::string count_line;
  int num_of_lines = 0;
  while( std::getline(linecounter, count_line ) ) num_of_lines += 1;

  cout<<sv_rcr_file<<" has "<<num_of_lines<<" lines in total.\n";
  linecounter.close();

  // We assume each rcr outlet has 6 lines.
  int num_of_chunk = num_of_lines - 1;
  
  if( num_of_chunk % 6 != 0 ) SYS_T::print_fatal("Error: %s number of lines wrong.\n", sv_rcr_file.c_str() );

  if( num_of_chunk / 6 != num_outlet ) SYS_T::print_fatal("Error: %s number of lines does not match num_outlet.\n", sv_rcr_file.c_str() );

  // Now start reading the info from the file.
  std::ifstream reader;
  reader.open( sv_rcr_file.c_str(), std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;
  
  // First line should be 2
  std::getline(reader, sline);
  sstrm.str( sline );
  int line_0;
  sstrm >> line_0;
  sstrm.clear();

  if( line_0 != 2 ) SYS_T::print_fatal( "Error: %s format is wrong. The first lineshould be 2.\n", sv_rcr_file.c_str() );

  std::vector<double> Rp, C, Rd, Pd;
  Rp.clear(); C.clear(); Rd.clear(); Pd.clear();

  for(int ii=0; ii<num_outlet; ++ii)
  {
    std::getline(reader, sline);
    
    if(reader.eof()) SYS_T::print_fatal("Error: EOF reached at line 1. Possibly the num_outlet is wrong.\n");

    sstrm.str( sline );
    int a;
    sstrm >> a;
    sstrm.clear();

    if( a != 2 ) SYS_T::print_fatal( "Error: numDataRCR has to be 2.\n" );

    // Rp
    std::getline( reader, sline );
    
    if(reader.eof()) SYS_T::print_fatal("Error: EOF reached at line 2. Possibly the num_outlet is wrong.\n");
    
    sstrm.str( sline );
    
    double temp_Rp;
    sstrm >> temp_Rp;
    sstrm.clear();

    Rp.push_back( temp_Rp );

    // C
    std::getline( reader, sline );
    
    if(reader.eof()) SYS_T::print_fatal("Error: EOF reached at line 3. Possibly the num_outlet is wrong.\n");
    
    sstrm.str( sline );

    double temp_C;
    sstrm >> temp_C;
    sstrm.clear();

    C.push_back( temp_C );

    // Rd
    std::getline( reader, sline );
    
    if(reader.eof()) SYS_T::print_fatal("Error: EOF reached at line 4. Possibly the num_outlet is wrong.\n");
    
    sstrm.str( sline );

    double temp_Rd;
    sstrm >> temp_Rd;
    sstrm.clear();

    Rd.push_back( temp_Rd );

    // Pd line 1
    std::getline( reader, sline );
    
    if(reader.eof()) SYS_T::print_fatal("Error: EOF reached at line 5. Possibly the num_outlet is wrong.\n");
    
    sstrm.str( sline );
    double temp_t1, temp_Pd;
    sstrm >> temp_t1;
    sstrm >> temp_Pd;
    sstrm.clear();

    Pd.push_back( temp_Pd );

    // Pd line 2
    std::getline( reader, sline );
    
    if(reader.eof()) SYS_T::print_fatal("Error: EOF reached at line 6. Possibly the num_outlet is wrong.\n");
    
    sstrm.str( sline );
    double temp_t2, temp_Pd_2;
    sstrm >> temp_t2;
    sstrm >> temp_Pd_2;
    sstrm.clear();

    if( temp_Pd != temp_Pd_2 ) SYS_T::print_fatal("Error: Pd_1 != Pd_2 for outlet %i \n", ii);
  }
  
  std::getline( reader, sline );

  if( !reader.eof() ) SYS_T::print_fatal("Error: file EOF is NOT reached. Possibly the num_outlet is samller than the number from the file.\n");

  reader.close();
  
  cout<<"Status: Finished reading "<<sv_rcr_file<<endl;

  // Start writing file
  std::ofstream ofile;
  ofile.open( out_rcr_file.c_str(), std::ofstream::out | std::ofstream::trunc );  
  
  ofile<<"RCR\t"<<num_outlet<<endl<<endl;

  for(int ii=0; ii<num_outlet; ++ii)
    ofile<<ii<<'\t'<<Rp[ii]<<'\t'<<C[ii]<<'\t'<<Rd[ii]<<'\t'<<Pd[ii]<<endl;

  ofile.close();

  cout<<"Status: Finished writing "<<out_rcr_file<<endl;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
