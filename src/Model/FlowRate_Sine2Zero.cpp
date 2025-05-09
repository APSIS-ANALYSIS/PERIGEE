#include "FlowRate_Sine2Zero.hpp"

FlowRate_Sine2Zero::FlowRate_Sine2Zero( const std::string &filename )
{
  SYS_T::commPrint("FlowRate_Sine2Zero: data read from %s \n", filename.c_str() );

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

  if( bc_type.compare("Sine") == 0 || bc_type.compare("SINE") == 0 || bc_type.compare("sine") == 0 )
  {
    coef_a.resize(num_nbc); coef_b.resize(num_nbc);
    num_of_mode.resize(num_nbc); w.resize(num_nbc); period.resize(num_nbc);
    thred_time.resize(num_nbc);
    start_flow_rate.resize(num_nbc);
    TI_std_dev.resize(num_nbc);
  }
  else
    SYS_T::print_fatal("FlowRate_Sine2Zero Error: inlet BC type in %s should be Inflow.\n", filename.c_str());

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

        if( face_id != nbc_id ) SYS_T::print_fatal("FlowRate_Sine2Zero Error: nbc in %s should be listed in ascending order.\n", filename.c_str());

        if(!(sstrm >> num_of_mode[nbc_id])) SYS_T::print_fatal("FlowRate_Sine2Zero Error: num_of_mode of nbc %d is undefined!\n", nbc_id);
        
        if(!(sstrm >> w[nbc_id])) SYS_T::print_fatal("FlowRate_Sine2Zero Error: w of nbc %d is undefined!\n", nbc_id);
        
        if(!(sstrm >> period[nbc_id])) SYS_T::print_fatal("FlowRate_Sine2Zero Error: period of nbc %d is undefined!\n", nbc_id);

        if(!(sstrm >> thred_time[nbc_id])) SYS_T::print_fatal("FlowRate_Sine2Zero Error: thred_time of nbc %d is undefined!\n", nbc_id);

        if(!(sstrm >> start_flow_rate[nbc_id]))
        {
          start_flow_rate[nbc_id] = 0.0;
          SYS_T::commPrint("FlowRate_Sine2Zero: default start_flow_rate of nbc %d is 0.0\n", nbc_id);
        }

        if(!(sstrm >> TI_std_dev[nbc_id]))
        {
          TI_std_dev[nbc_id] = 0.0;
          SYS_T::commPrint("FlowRate_Sine2Zero: default TI_std_dev of nbc %d is 0.0\n", nbc_id);
        }

        sstrm.clear();
        break;
      }
    }

    // Check the compatibility of period and w. If the difference
    // is larger than 0.01, print a warning message
    if( std::abs(2.0 * MATH_T::PI / period[nbc_id] - w[nbc_id] ) >= 0.01 ) SYS_T::commPrint( "\nFlowRate_Sine2Zero WARNING: nbc_id %d incompatible period and w, \n2xpi/period = %e and w = %e.\n", nbc_id, 2.0*static_cast<double>(MATH_T::PI)/period[nbc_id], w[nbc_id] );

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

    coef_a[nbc_id].shrink_to_fit();
    
    if( static_cast<int>(coef_a[nbc_id].size()) != num_of_mode[nbc_id]+1 )
      SYS_T::print_fatal("FlowRate_Sine2Zero Error: nbc_id %d a-coefficients in %s incompatible with the given number of modes.\n", nbc_id, filename.c_str());

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

    coef_b[nbc_id].shrink_to_fit();

    if( static_cast<int>(coef_b[nbc_id].size()) != num_of_mode[nbc_id]+1 )
      SYS_T::print_fatal("FlowRate_Sine2Zero Error: nbc_id %d b-coefficients in %s incompatible with the given number of modes.\n", nbc_id, filename.c_str());
  }

  // Finish reading the file and close it
  reader.close();

  // Generate the target flow rate as the sum of coef_a.
  target_flow_rate.resize( num_nbc );

  for(int ii=0; ii<num_nbc; ++ii) target_flow_rate[ii] = VEC_T::sum( coef_a[ii] );

  MPI_Barrier(PETSC_COMM_WORLD);
}

double FlowRate_Sine2Zero::get_flow_rate( const int &nbc_id,
    const double &time ) const
{
  double out_dot_rate = 0.0;

  if( time < thred_time[nbc_id] && time >= 0.0 ) 
    out_dot_rate = 0.5 * (target_flow_rate[nbc_id]- start_flow_rate[nbc_id]) * MATH_T::PI / thred_time[nbc_id] * std::sin(MATH_T::PI * time / thred_time[nbc_id]);

  return out_dot_rate;
}

void FlowRate_Sine2Zero::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("     FlowRate_Sine2Zero:\n");
  for(int nbc_id = 0; nbc_id < num_nbc; ++nbc_id)
  {
    SYS_T::commPrint("  -- nbc_id = %d \n", nbc_id);
    SYS_T::commPrint("     start flow rate : %e \n", start_flow_rate[nbc_id]);
    SYS_T::commPrint("     time to reach steady state : %e \n", thred_time[nbc_id]);
    SYS_T::commPrint("     target flow rate : %e \n", target_flow_rate[nbc_id]);
    SYS_T::commPrint("     turbulance intensity : %e \n", TI_std_dev[nbc_id]);
  }
  SYS_T::print_sep_line();
}

// EOF
