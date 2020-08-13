#include "Sys_Tools.hpp"

void SYS_T::print_perigee_art()
{
  PetscPrintf(PETSC_COMM_WORLD, "$$$$$$$\\  $$$$$$$$\\ $$$$$$$\\  $$$$$$\\  $$$$$$\\  $$$$$$$$\\ $$$$$$$$\\ \n");
  PetscPrintf(PETSC_COMM_WORLD, "$$  __$$\\ $$  _____|$$  __$$\\ \\_$$  _|$$  __$$\\ $$  _____|$$  _____| \n");
  PetscPrintf(PETSC_COMM_WORLD, "$$ |  $$ |$$ |      $$ |  $$ |  $$ |  $$ /  \\__|$$ |      $$ | \n");
  PetscPrintf(PETSC_COMM_WORLD, "$$$$$$$  |$$$$$\\    $$$$$$$  |  $$ |  $$ |$$$$\\ $$$$$\\    $$$$$\\ \n");
  PetscPrintf(PETSC_COMM_WORLD, "$$  ____/ $$  __|   $$  __$$<   $$ |  $$ |\\_$$ |$$  __|   $$  __| \n");
  PetscPrintf(PETSC_COMM_WORLD, "$$ |      $$ |      $$ |  $$ |  $$ |  $$ |  $$ |$$ |      $$ | \n");
  PetscPrintf(PETSC_COMM_WORLD, "$$ |      $$$$$$$$\\ $$ |  $$ |$$$$$$\\ \\$$$$$$  |$$$$$$$$\\ $$$$$$$$\\ \n");
  PetscPrintf(PETSC_COMM_WORLD, "\\__|      \\________|\\__|  \\__|\\______| \\______/ \\________|\\________| \n \n");
}


int SYS_T::get_genbc_file_type( const char * const &lpn_filename )
{
  const std::string temp_name( lpn_filename );
  SYS_T::file_check( temp_name );

  std::ifstream reader;
  reader.open( lpn_filename, std::ifstream::in );

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
  else
  {
    return 0;
  }
}


double SYS_T::gen_randomD_closed(const double &min, const double &max)
{
  return ( rand() % 1000001 ) * 1.0e-6 * (max - min) + min;
}


double SYS_T::gen_randomD_open(const double &min, const double &max)
{
  return ( rand() % 999998 + 1 ) * 1.0e-6 * (max - min) + min;
}

int SYS_T::gen_randomI_closed(const int &min, const int &max)
{
  return ( rand() % (max - min + 1)) + min;
}


void SYS_T::print_MaxMemUsage()
{
  PetscLogDouble memo = 0.0, memototal = 0.0;
  PetscMemoryGetMaximumUsage(&memo);
  MPI_Reduce(&memo, &memototal, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(rank == 0)
  {
    std::cout<<"\n Maximum Memeory usage : ";
    print_mem_size(memototal); std::cout<<"\n";
  }
}


void SYS_T::print_CurMemUsage()
{
  PetscLogDouble memo = 0.0, memototal = 0.0;
  PetscMemoryGetCurrentUsage(&memo);
  MPI_Reduce(&memo, &memototal, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(rank == 0)
  {
    std::cout<<"\n Current Memeory usage : ";
    print_mem_size(memototal); std::cout<<"\n";
  }
}


void SYS_T::print_MaxMallocUsage()
{
  PetscLogDouble memo = 0.0, memototal = 0.0;
  PetscMallocGetMaximumUsage(&memo);
  MPI_Reduce(&memo, &memototal, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(rank == 0)
  {
    std::cout<<"\n Maximum PETSc malloced : ";
    print_mem_size(memototal); std::cout<<"\n";
  }
}


void SYS_T::print_CurMallocUsage()
{
  PetscLogDouble memo = 0.0, memototal = 0.0;
  PetscMallocGetCurrentUsage(&memo);
  MPI_Reduce(&memo, &memototal, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(rank == 0)
  {
    std::cout<<"\n Current PETSc malloced : ";
    print_mem_size(memototal); std::cout<<"\n";
  }
}


std::string SYS_T::get_time()
{
  std::time_t  time1= std::time (0);
  std::tm     *time = std::localtime(&time1);

  std::ostringstream o;
  o << time->tm_hour << ":"
    << (time->tm_min < 10 ? "0" : "") << time->tm_min << ":"
    << (time->tm_sec < 10 ? "0" : "") << time->tm_sec;

  return o.str();
}


std::string SYS_T::get_date()
{
  std::time_t  time1= std::time (0);
  std::tm     *time = std::localtime(&time1);

  std::ostringstream o;
  o << time->tm_year + 1900 << "/"
    << time->tm_mon + 1 << "/"
    << time->tm_mday;

  return o.str();
}


void SYS_T::get_memory_stats (MemoryStats &stats)
{
  stats.VmPeak = stats.VmSize = stats.VmHWM = stats.VmRSS = 0;

#if defined(__linux__)
  std::ifstream file("/proc/self/status");
  std::string line;
  std::string name;
  while (!file.eof())
  {
    file >> name;
    if (name == "VmPeak:")
      file >> stats.VmPeak;
    else if (name == "VmSize:")
      file >> stats.VmSize;
    else if (name == "VmHWM:")
      file >> stats.VmHWM;
    else if (name == "VmRSS:")
    {
      file >> stats.VmRSS;
      break; //this is always the last entry
    }

    getline(file, line);
  }
#endif
}


SYS_T::Timer::Timer()
{
  startedAt = 0;
  stoppedAt = 0;
}

SYS_T::Timer::~Timer()
{}

void SYS_T::Timer::Start()
{
  startedAt = clock();
}


void SYS_T::Timer::Stop()
{
  stoppedAt = clock();
}


void SYS_T::Timer::Reset()
{
  startedAt = 0;
  stoppedAt = 0;
}


double SYS_T::Timer::get_sec() const
{
  return (double)(stoppedAt - startedAt)/(double)CLOCKS_PER_SEC;
}

// EOF
