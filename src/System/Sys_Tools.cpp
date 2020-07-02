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


int SYS_T::readFile(std::ifstream& infile,
    std::vector<double>& sKnots,
    std::vector<double>& tKnots,
    std::vector<double>& uKnots,
    int& sDegree, int& tDegree,
    int& uDegree, int& numCPts,
    std::vector<double>& ctrlPts)
{
  std::istringstream sstrm;
  std::string line;
  std::string sWord;
  char inWord[1001];
  char x;

  int hits = 0;
  while( hits < 8 && !infile.eof() )
  {
    infile.getline(inWord, 1000);
    line = inWord;
    sstrm.str( line );
    sstrm >> inWord;

    sWord = inWord;
    if( sWord == "TYPE" )
    {
      sstrm >> x;
      if( x != '=' )
      {
        std::cerr << "Bad file format: 1" << std::endl;
        return 1;
      }
      sstrm >> inWord;
      sWord = inWord;
      if( sWord != "NURBS" )
      {
        std::cerr << "The file format is not NURBS.  Cannot read file." << std::endl;
        return 1;
      }

      hits++;
    }
    else if( sWord == "GLOBAL_S" )
    {
      SYS_T::readKnotVector(sstrm, sKnots);
      hits++;
    }
    else if( sWord == "GLOBAL_T" )
    {
      SYS_T::readKnotVector(sstrm, tKnots);
      hits++;
    }    
    else if( sWord == "GLOBAL_U" )
    {
      SYS_T::readKnotVector(sstrm, uKnots);
      hits++;
    }
    else if( sWord == "DEGREE_S" )
    {
      sstrm >> x;
      if( x != '=' || sstrm.eof() )
      {
        std::cerr << "Bad file format: 3.1" << std::endl;
        return 1;
      }
      sstrm >> sDegree;
      hits++;
    }
    else if( sWord == "DEGREE_T" )
    {
      sstrm >> x;
      if( x != '=' || sstrm.eof() )
      {
        std::cerr << "Bad file format: 3.2" << std::endl;
        return 1;
      }
      sstrm >> tDegree;
      hits++;
    }    
    else if( sWord == "DEGREE_U" )
    {
      sstrm >> x;
      if( x != '=' || sstrm.eof() )
      {
        std::cerr << "Bad file format: 3.3" << std::endl;
        return 1;
      }
      sstrm >> uDegree;
      hits++;
    }
    else if( sWord == "NUM_CP" )
    {
      sstrm >> x;
      if( x != '=' || sstrm.eof() )
      {
        std::cerr << "Bad file format: 4.1" << std::endl;
        return 1;
      }
      sstrm >> numCPts;
      hits++;
    }

    sstrm.clear();
  }
  if( 8 != hits )
  {
    std::cerr << "Error reading file.  Something is wrong with the format\n" << std::endl;
    return 1;
  }

  int ii;
  for( ii = 0; ii < numCPts*4; ii++ )
  {
    if( infile.eof() )
    {
      std::cerr << "Wrong number of control points." << std::endl;
      return 1;
    }

    double coord;
    infile >> coord;
    ctrlPts.push_back( coord );
  }

  if( (int) ctrlPts.size() != numCPts*4 )
  {
    std::cerr << "Error reading control points" << std::endl
      << ctrlPts.size() << std::endl;
    return 1;
  }

  return 0;
};



int SYS_T::readSHFile( std::ifstream &infile,
    std::vector<double> &sKnots,
    std::vector<double> &tKnots,
    std::vector<double> &uKnots,
    int &sDegree, int &tDegree,
    int &uDegree, int &numCPts,
    std::vector<double> &ctrlPts )
{
  std::istringstream sstrm;
  std::string sline;

  getline(infile, sline); // 3
  sstrm.clear();

  getline(infile, sline);
  sstrm.str(sline);
  sstrm>>sDegree;
  sstrm>>tDegree;
  sstrm>>uDegree;
  sstrm.clear();

  getline(infile, sline);
  sstrm.str(sline);
  int len_s, len_t, len_u;
  sstrm>>len_s; 
  sstrm>>len_t;
  sstrm>>len_u;
  sstrm.clear();

  getline(infile, sline);
  sstrm.str(sline);
  double knot;
  while(!sstrm.eof())
  {
    sstrm>>knot;
    sKnots.push_back( knot );
  }
  sstrm.clear();

  getline(infile, sline);
  sstrm.str(sline);
  while(!sstrm.eof())
  {
    sstrm>>knot;
    tKnots.push_back( knot );
  }
  sstrm.clear();

  getline(infile, sline);
  sstrm.str(sline);
  while(!sstrm.eof())
  {
    sstrm>>knot;
    uKnots.push_back( knot );
  }
  sstrm.clear();

  // Now we check if the given number of basis function is compatible with the
  // given knot vector
  if( (int) sKnots.size() != len_s + sDegree + 1 )
  {
    std::cerr<<" Error: sKnot size is incompatible with the number of basis in s-direction. \n";
    return 1;
  }

  if( (int) tKnots.size() != len_t + tDegree + 1 )
  {
    std::cerr<<" Error: tKnot size is incompatible with the number of basis in t-direction. \n";
    return 1;
  }

  if( (int) uKnots.size() != len_u + uDegree + 1 )
  {
    std::cerr<<" Error: uKnot size is incompatible with the number of basis in u-direction. \n";
    return 1;
  }
  
  numCPts = len_s * len_t * len_u;
  double in_cp;
  for(int ii=0; ii<numCPts; ++ii)
  {
    getline(infile, sline); sstrm.str(sline);
    sstrm>>in_cp; ctrlPts.push_back(in_cp);
    sstrm>>in_cp; ctrlPts.push_back(in_cp);
    sstrm>>in_cp; ctrlPts.push_back(in_cp);
    sstrm>>in_cp; ctrlPts.push_back(in_cp);
    
    if(!sstrm.eof())
    {
      std::cerr<<"Error: wrong format in control points. \n";
      return 1;
    }
    
    sstrm.clear();
  }

  return 0;
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


void SYS_T::checkSHKnots( std::vector<double> &sKnots, const int &sdeg )
{
  const int len = sKnots.size();

  const double fr = sKnots[0];
  const double ba = sKnots[len-1];

  for(int ii=0; ii<sdeg+1; ++ii)
  {
    if( sKnots[ii] != fr || sKnots[len-1-ii] != ba )
    {
      std::cerr<<"Error: knot vector is incompatible with the poly degree. ";
      std::cerr<<"There is not enough knot to make open knot vector. \n";
      exit(1);
    }
  }

  if( sKnots[sdeg+1] == fr || sKnots[len-sdeg-2] == ba )
  {
    std::cerr<<"Error: knot vector is incompatible with the poly degree ";
    std::cerr<<" Extra repeated knot found. \n";
    exit(1);
  }

  if( fr != 0.0 || ba != 1.0 )
  {
    const double invh = 1.0 / ( ba - fr );
    for(int ii=0; ii<len; ++ii)
      sKnots[ii] = (sKnots[ii] - fr) * invh ;
  }

  for(int ii=0; ii<sdeg+1; ++ii)
  {
    sKnots[ii] = 0.0;
    sKnots[len-1-ii] = 1.0;
  }
}



int SYS_T::readKnotVector(std::istringstream& sstrm,
    std::vector<double>& knots)
{
  char x;

  //Discard = and [ characters
  sstrm >> x;
  if( x != '=' )
  {
    std::cerr << "Bad file format: 2.1" << std::endl;
    return 1;
  }
  sstrm >> x;
  if( x != '[' )
  {
    std::cerr << "Bad file format: 2.2" << std::endl;
    return 1;
  }

  sstrm >> x;
  while( x != ']' && !sstrm.eof() )
  {
    sstrm.putback( x );
    double knot;
    sstrm >> knot;
    knots.push_back( knot );
    sstrm >> x;
  }

  return 0;
};


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
