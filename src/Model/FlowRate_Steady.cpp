#include "FlowRate_Steady.hpp"

FlowRate_Steady::FlowRate_Steady( const std::string &filename )
{
  SYS_T::commPrint("FlowRate_Steady: data read from %s \n", filename.c_str() );

  SYS_T::file_check( filename );
  
  std::ifstream reader;
  reader.open( filename.c_str(), std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;
  std::string bc_type;

  // The first non-commented, non-empty line should be
  // Inflow num_nbc
  while( std::getline(reader, sline) )
  {
    if( sline[0] !='#' && !sline.empty() )
    {
      sstrm.str(sline);
      sstrm >> bc_type;
      sstrm >> num_nbc;
      sstrm.clear();
      break;
    }
  }

  std::vector< std::vector<double> > coef_a, coef_b;
  std::vector<int> num_of_mode;
  std::vector<double> w, period;

  if( bc_type.compare("Steady") == 0 || bc_type.compare("STEADY") == 0 || bc_type.compare("steady") == 0 )
  {
    coef_a.resize(num_nbc); coef_b.resize(num_nbc);
    num_of_mode.resize(num_nbc); w.resize(num_nbc); period.resize(num_nbc);
  }
  else
    SYS_T::print_fatal("FlowRate_Steady Error: inlet BC type in %s should be Inflow.\n", filename.c_str());

  // Read in num_of_mode, w, period, coef_a, and coef_b per nbc
  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    while( std::getline(reader, sline) )
    {
      // face_id num_of_mode w period
      if( sline[0] !='#' && !sline.empty() )
      {
        sstrm.str(sline);
        int face_id;
        sstrm >> face_id;

        if( face_id != nbc_id ) SYS_T::print_fatal("FlowRate_Steady Error: nbc in %s should be listed in ascending order.\n", filename.c_str());

        sstrm >> num_of_mode[nbc_id];
        sstrm >> w[nbc_id];
        sstrm >> period[nbc_id];

        sstrm.clear();
        break;
      }
    }

    // Check the compatibility of period and w. If the difference
    // is larger than 0.01, print a warning message
    if( std::abs(2.0 * MATH_T::PI / period[nbc_id] - w[nbc_id] ) >= 0.01 )
      SYS_T::commPrint( "\nFlowRate_Steady WARNING: nbc_id %d incompatible period and w, \n2xpi/period = %e and w = %e.\n", nbc_id, 2.0*MATH_T::PI/period[nbc_id], w[nbc_id] );

    coef_a[nbc_id].clear(); coef_b[nbc_id].clear(); 

    while( std::getline(reader, sline) )
    {
      // coef_a
      if( sline[0] !='#' && !sline.empty() )
      {
        sstrm.str(sline);
        double temp_coef;
        while( sstrm >> temp_coef ) coef_a[nbc_id].push_back( temp_coef );
     
        sstrm.clear(); 
        break;
      }
    }

    VEC_T::shrink2fit( coef_a[nbc_id] );
    
    if( static_cast<int>(coef_a[nbc_id].size()) != num_of_mode[nbc_id]+1 )
      SYS_T::print_fatal("FlowRate_Steady Error: nbc_id %d a-coefficients in %s incompatible with the given number of modes.\n", nbc_id, filename.c_str());

    while( std::getline(reader, sline) )
    {
      // coef_b
      if( sline[0] !='#' && !sline.empty() )
      {
        sstrm.str(sline);
        double temp_coef;
        while( sstrm >> temp_coef ) coef_b[nbc_id].push_back( temp_coef );
     
        sstrm.clear(); 
        break;
      }
    }

    VEC_T::shrink2fit( coef_b[nbc_id] );

    if( static_cast<int>(coef_b[nbc_id].size()) != num_of_mode[nbc_id]+1 )
      SYS_T::print_fatal("FlowRate_Steady Error: nbc_id %d b-coefficients in %s incompatible with the given number of modes.\n", nbc_id, filename.c_str());
  }

  // Finish reading the file and close it
  reader.close();

  // Define flow rate to be the smallest value in flow_waveform
  flowrate.resize( num_nbc );

  for(int ii=0; ii<num_nbc; ++ii) flowrate[ii] = VEC_T::sum( coef_a[ii] );

  // Calculate the flow rate and record them on disk as 
  // Inlet_XXX_flowrate.txt with sampling interval 0.001
  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    if( SYS_T::get_MPI_rank() == 0 )
    {
      std::ofstream ofile;
      ofile.open( gen_flowfile_name(nbc_id).c_str(), std::ofstream::out | std::ofstream::trunc );
      const double tt = 0.0;
      ofile << tt <<'\t'<<get_flow_rate(nbc_id, tt)<< '\n';
      ofile.close();
    }
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

double FlowRate_Steady::get_flow_rate(const int &nbc_id , const double &time) const
{
  return flowrate[nbc_id];
}

void FlowRate_Steady::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("  FlowRate_Steady:\n");
  for(int nbc_id = 0; nbc_id < num_nbc; ++nbc_id)
  {
    SYS_T::commPrint("  -- nbc_id = %d", nbc_id);
    SYS_T::commPrint("     flow rate =%e \n", flowrate[nbc_id]);
  }
  SYS_T::print_sep_line();
}

// EOF
