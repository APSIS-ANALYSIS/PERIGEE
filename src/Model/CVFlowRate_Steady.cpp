#include "CVFlowRate_Steady.hpp"

CVFlowRate_Steady::CVFlowRate_Steady( const char * const &filename )
{
  SYS_T::commPrint("CVFlowRate_Steady: data read from %s \n", filename);

  const std::string temp_name(filename);
  SYS_T::file_check( temp_name );

  std::ifstream reader;
  reader.open( filename, std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;

  int num_of_mode;
  double w, period;
  while( std::getline(reader, sline) )
  {
    // The first non-commented non-empty line records num-of-mode & w & period
    if( sline[0] !='#' && !sline.empty() )
    {
      sstrm.str(sline);
      sstrm >> num_of_mode;
      sstrm >> w;
      sstrm >> period;
      sstrm.clear();
      break;
    }
  }

  std::vector<double> coef_a, coef_b;
  // Allocate coef_a & coef_b
  coef_a.clear(); coef_b.clear();

  while( std::getline(reader, sline) )
  {
    // The second non-commented non-empty line records coef_a
    if( sline[0] !='#' && !sline.empty() )
    {
      sstrm.str(sline);
      double temp_coef;
      while( sstrm >> temp_coef ) coef_a.push_back( temp_coef );

      sstrm.clear();
      break;
    }
  }

  VEC_T::shrink2fit(coef_a);

  if(static_cast<int>(coef_a.size()) != num_of_mode+1) SYS_T::print_fatal("Error: CVFlowRate_Steady input file %s, the a-coefficient array is incompatible with the given number of modes.\n", filename);

  while( std::getline(reader, sline) )
  {
    // The third non-commented non-empty line records coef_b
    if( sline[0] !='#' && !sline.empty() )
    {
      sstrm.str(sline);
      double temp_coef;
      while( sstrm >> temp_coef ) coef_b.push_back( temp_coef );

      sstrm.clear();
      break;
    }
  }

  VEC_T::shrink2fit(coef_b);

  if(static_cast<int>(coef_b.size()) != num_of_mode+1) SYS_T::print_fatal("Error: CVFlowRate_Steady input file %s, the b-coefficient array is incompatible with the given number of modes.\n", filename);

  // Finish reading the file and close it
  reader.close();

  // Check the compatibility of period and w, if the difference
  // is larger than 0.01, print a warning message
  if( std::abs(2.0 * MATH_T::PI / period - w ) >= 0.01 ) SYS_T::commPrint("\nWARNING: CVFlowRate_Steady period and w does not match well: \n2xpi/period = %e and w = %e.\n", 2.0*MATH_T::PI/period, w);

  std::vector<double> flow_waveform; flow_waveform.resize( static_cast<int>(period / 0.001) + 1 );
  for(unsigned int ii = 0; ii < flow_waveform.size(); ++ii )
  {
    flow_waveform[ii] = coef_a[0];
    
    for(int jj=1; jj<=num_of_mode; ++jj)
      flow_waveform[ii] += coef_a[jj]*cos(jj*w*ii*0.001) + coef_b[jj]*sin(jj*w*ii*0.001);
  }

  // Define flow rate to be the smallest value in flow_waveform
  flowrate = *std::min_element(flow_waveform.begin(), flow_waveform.end());
  
  // Calculate the flow rate and record them on disk as
  // Inlet_flowrate.txt with sampling interval 0.001
  if( SYS_T::get_MPI_rank() == 0 )
  {
    std::ofstream ofile;
    ofile.open( "Inlet_flowrate.txt", std::ofstream::out | std::ofstream::trunc );
    for(unsigned int ii = 0; ii < flow_waveform.size(); ++ii )
      ofile << ii * 0.001 << "\t" << flowrate << "\n";
    ofile.close();
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

CVFlowRate_Steady::CVFlowRate_Steady( const double &in_flowrate ) : flowrate( in_flowrate )
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    std::ofstream ofile;
    ofile.open( "Inlet_flowrate.txt", std::ofstream::out | std::ofstream::trunc );
    ofile <<"constant flow rate is "<< flowrate << "\n";
    ofile.close();
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

CVFlowRate_Steady::~CVFlowRate_Steady()
{
}

double CVFlowRate_Steady::get_flow_rate(const double &time) const
{
  return flowrate;
}

void CVFlowRate_Steady::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("  CVFlowRate_Steady:\n");
  SYS_T::commPrint("  flowrate = %e \n", flowrate);
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

// EOF
