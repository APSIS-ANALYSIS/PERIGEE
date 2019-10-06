#include "CVFlowRate_Unsteady.hpp"

CVFlowRate_Unsteady::CVFlowRate_Unsteady( const char * const &filename )
{
  SYS_T::commPrint("CVFlowRate_Unsteady: data read from %s \n", filename);

  const std::string temp_name(filename);
  SYS_T::file_check( temp_name );
  
  std::ifstream reader;
  reader.open( filename, std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;

  while( std::getline(reader, sline) )
  {
    // The first non-commented non-empty line records
    // num-of-mode w period
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
  
  if(static_cast<int>(coef_a.size()) != num_of_mode+1) SYS_T::print_fatal("Error: CVFlowRate_Unsteady input file %s, the a-coefficient array is incompatible with the given number of modes.\n", filename);

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

  if(static_cast<int>(coef_b.size()) != num_of_mode+1) SYS_T::print_fatal("Error: CVFlowRate_Unsteady input file %s, the b-coefficient array is incompatible with the given number of modes.\n", filename);

  // Finish reading the file and close it
  reader.close();

  // Check the compatibility of period and w, if the difference
  // is larger than 0.01, print a warning message
  if( std::abs(2.0 * MATH_T::PI / period - w ) >= 0.01 ) SYS_T::commPrint("Warning: CVFlowRate_Unsteady period and w does not match well: \n2xpi/period = %e and w = %e.\n", 2.0*MATH_T::PI/period, w);

  // Calculate the flow rate and record them on disk as 
  // Inlet_flowrate.txt with sampling interval 0.001
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if( rank == 0 )
  {
    std::ofstream ofile;
    ofile.open( "Inlet_flowrate.txt", std::ofstream::out | std::ofstream::trunc );
    for(double tt=0; tt <= period; tt += 0.001 )
      ofile<<tt<<'\t'<<get_flow_rate(tt)<<'\n';
    ofile.close();
  }
}


CVFlowRate_Unsteady::~CVFlowRate_Unsteady()
{
  VEC_T::clean(coef_a); VEC_T::clean(coef_b);
}


double CVFlowRate_Unsteady::get_flow_rate(const double &time) const
{
  const int num_of_past_period = time / period;
  const double local_time = time - num_of_past_period * period;

  double sum = coef_a[0];
  for(int ii=1; ii<=num_of_mode; ++ii)
    sum += coef_a[ii]*cos( ii*w*local_time ) + coef_b[ii]*sin(ii*w*local_time);

  return sum;
}

void CVFlowRate_Unsteady::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("     CVFlowRate_Unsteady: ");
  PetscPrintf(PETSC_COMM_WORLD, " w = %e, period =%e \n", w, period);
  PetscPrintf(PETSC_COMM_WORLD, "     a[0] + Sum{ a[i] cos(i x w x t) + b[i] sin(i x w x t) }, for i = 1,...,%d. \n", num_of_mode);
  for(int ii=0; ii<=num_of_mode; ++ii) PetscPrintf(PETSC_COMM_WORLD, "     i = %d, a = %e, b = %e \n", ii, coef_a[ii], coef_b[ii]);
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

// EOF
