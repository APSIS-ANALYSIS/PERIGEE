#ifndef SYS_TOOLS_HPP
#define SYS_TOOLS_HPP
// ==================================================================
// Sys_Tools.hpp
// ------------------------------------------------------------------
// The SYS_T namespace contains a suite of tools at the system level.
//
// These functions will be frequently used in the PERIGEE code.
// ==================================================================
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include "petsc.h"

namespace SYS_T
{
  // Print ASCII art fonts
  void print_perigee_art();

  // Return the rank of the CPU
  inline PetscMPIInt get_MPI_rank()
  {
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    return rank;
  }

  // Return the number of total number of CPUs
  inline PetscMPIInt get_MPI_size()
  {
    PetscMPIInt size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    return size;
  }

  // ----------------------------------------------------------------
  // for Zinan, replace this by GenBC_T::get_genbc_file_type
  // ! get_genbc_file_type : read the genbc file and determine what
  //   type of gen bc the file specifies. It will return   
  //   0 for unknown type
  //   1 for Resistance
  //   2 for RCR
  // ----------------------------------------------------------------
  int get_genbc_file_type( const char * const &lpn_filename );

  // ----------------------------------------------------------------
  // gen_partfile_name( baseName, rank )
  // Generate a partition file's name (hdf5 file) in the default
  // manner. It will return baseName_pxxxxx.h5.
  // ----------------------------------------------------------------
  inline std::string gen_partfile_name( const std::string &baseName, 
      const int &rank )
  {
    std::ostringstream ss;
    ss << baseName <<"_p";

    if( rank / 10 == 0 )
      ss<<"0000";
    else if( rank / 100 == 0 )
      ss<<"000";
    else if( rank / 1000 == 0 )
      ss<<"00";
    else if( rank / 10000 == 0 )
      ss<<"0";

    ss<<rank<<".h5";
    return ss.str();
  }

  // ----------------------------------------------------------------
  // get_xyz_index() 
  // Assume ii = iz * dim_x * dim_y + iy * dim_x + ix
  // this function will return ix iy and iz based on the input ii, 
  // dim_x, dim_y.
  // ----------------------------------------------------------------
  inline void get_xyz_index( const int &ii, const int &dim_x,
      const int &dim_y, int &ix, int &iy, int &iz)
  {
    const int ixy = ii % (dim_x * dim_y);
    iz = (ii - ixy) / (dim_x * dim_y);
    ix = ixy % dim_x; iy = (ixy - ix) / dim_x;
  }

  // ----------------------------------------------------------------
  // get_xy_index()
  // Assume ii = iy * dim_x + ix;
  // this function will return ix and iy based on the input ii and dim_x.
  // ----------------------------------------------------------------
  inline void get_xy_index( const int &ii, const int &dim_x, int &ix, int &iy)
  {
    ix = ii % dim_x;
    iy = (ii-ix)/dim_x;
  }

  // ----------------------------------------------------------------
  // gen_random()
  // Generate a random double in [min, max] domain for _closed; 
  // (min, max) for open.
  // Gernerate a random int in [min, max] domain.
  // NOTE: Users have to call srand(time(NULL)) before calling the 
  //       following three gen functions.
  // E.G.: srand(time(NULL));
  //       for-loop
  //       {gen_randomD_xxx(...); ...}
  // ----------------------------------------------------------------
  double gen_randomD_closed( const double &min, const double &max );

  double gen_randomD_open( const double &min, const double &max );

  int gen_randomI_closed( const int &min, const int &max );

  // ----------------------------------------------------------------
  // print size based on bytes
  // ----------------------------------------------------------------
  template<typename T> void print_mem_size( const T &byte_size )
  {
    if(byte_size > 1.0e9)
      std::cout<<double(byte_size)/1.0e9<<" GB.";
    else if(byte_size > 1.0e6)
      std::cout<<double(byte_size)/1.0e6<<" MB.";
    else if(byte_size > 1.0e3)
      std::cout<<double(byte_size)/1.0e3<<" KB";
    else
      std::cout<<byte_size<<" Bytes.";
  }

  template<typename T> std::string get_string_mem_size( const T &byte_size )
  {
    std::ostringstream ss;
    if(byte_size > 1.0e9)
      ss<<(double)byte_size/1.0e9<<" GB.";
    else if(byte_size > 1.0e6)
      ss<<(double)byte_size/1.0e6<<" MB.";
    else if(byte_size > 1.0e3)
      ss<<(double)byte_size/1.0e3<<" KB.";
    else
      ss<<(double)byte_size<<" B.";
    return ss.str();
  }

  // ----------------------------------------------------------------
  // to_string functions : convert numeric values to string 
  // Note: std::to_string is implemented in string in C++ 11.
  // ----------------------------------------------------------------
  inline std::string to_string( const int &a )
  {
    std::ostringstream ss;
    ss<<a;
    return ss.str();
  }

  inline std::string to_string( const double &a )
  {
    std::ostringstream ss;
    ss<<a;
    return ss.str();
  }

  // to get a non-constant (i.e. writable) char array from a string, 
  // use &b[0]. Note: string.c_str() returns a const char array.
  inline void to_char( const std::string &a, std::vector<char> &b )
  {
    b.assign(a.begin(), a.end());
    b.push_back('\0');
  }

  // this gets a non-constant/writable char array from a string.
  // this directly returns a char array, but the user is responsible 
  // for deleting the char pointer after usage.
  inline void to_char( const std::string &a, char * &b )
  {
    b = new char [a.size() + 1];
    std::copy(a.begin(), a.end(), b);
    b[a.size()] = '\0';
  }

  // ================================================================
  // The following are used in processor.
  // ================================================================
  // 1. Synchronized print on screen. Here we use petsc functions:
  //    PetscSynchronizedPrintf() and PetscSynchronizedFlush().
  //    The output should be lised from proc 0, 1, ... to n in sequence.
  inline void synPrint(const std::string &output, const int &cpu_rank)
  {
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Proc %d: ", cpu_rank);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, output.c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  }
  inline void synPrint(const char * const &output, const int &cpu_rank)
  {
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Proc %d: ", cpu_rank);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, output);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  }
  inline void synPrint(const std::string &output)
  {
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, output.c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  }
  inline void synPrint(const char * const &output)
  {
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, output);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  }

  // 2. print from processor 0, other preprocessors are ignored.
  //    PetscPrintf() with PETSC_COMM_WORLD is used.
  inline void commPrint(const char output[], ...)
  {
    PetscMPIInt rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if(!rank)
    {
      va_list Argp;
      va_start(Argp, output);
      (*PetscVFPrintf)(PETSC_STDOUT,output,Argp);
      va_end(Argp);
    }
  }

  // 3. Print the data name and its value on screen using PetscPrintf
  //    and PETSC_COMM_WORLD. This is particularly designed to print
  //    command line arguments.
  inline void cmdPrint(const char * const &dataname, const int &datavalue)
  {std::ostringstream ss; ss<<dataname<<" "<<datavalue<<"\n"; PetscPrintf(PETSC_COMM_WORLD, ss.str().c_str());}
  inline void cmdPrint(const char * const &dataname, const double &datavalue)
  {std::ostringstream ss; ss<<dataname<<" "<<datavalue<<"\n"; PetscPrintf(PETSC_COMM_WORLD, ss.str().c_str());}
  inline void cmdPrint(const char * const &dataname, const std::string &datavalue)
  {std::ostringstream ss; ss<<dataname<<" "<<datavalue<<"\n"; PetscPrintf(PETSC_COMM_WORLD, ss.str().c_str());}

  // 4. Print specific system message
  inline void synPrintElementInfo(int nlocalele, double memspace, 
      double totaltime, PetscMPIInt rank)
  {
    std::string memusage = get_string_mem_size(memspace);
    std::ostringstream ss; ss<<nlocalele<<" elements cached, "<<totaltime<<
      " secs and "<<memusage<<std::endl;
    synPrint(ss.str(), rank);
  }

  // 5. Print fatal error message and terminate the MPI process
  inline void print_fatal( const char output[], ... )
  {
    const PetscMPIInt rank = get_MPI_rank();

    if(!rank)
    {
      va_list Argp;
      va_start(Argp, output);
      (*PetscVFPrintf)(PETSC_STDOUT,output,Argp);
      va_end(Argp);
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  inline void print_fatal_if( bool a, const char output[], ... )
  {
    if( a )
    {
      const PetscMPIInt rank = get_MPI_rank();

      if(!rank)
      {
        va_list Argp;
        va_start(Argp, output);
        (*PetscVFPrintf)(PETSC_STDOUT,output,Argp);
        va_end(Argp);
      }

      MPI_Barrier(PETSC_COMM_WORLD);
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
  }

  // exit message printers are used in terminating serial program when 
  // the communicator for MPI is not available.
  inline void print_exit( const char * const &mesg )
  {
    std::cout<<mesg<<std::endl;
    exit( EXIT_FAILURE );
  }

  inline void print_exit( const std::string &mesg )
  {
    std::cout<<mesg<<std::endl;
    exit( EXIT_FAILURE );
  }

  inline void print_exit_if( bool a, const char * const &mesg )
  {
    if( a ) print_exit(mesg);
  }

  inline void print_exit_if( bool a, const std::string &mesg )
  {
    if( a ) print_exit(mesg);
  }

  // =================================================================
  // The followings are system function to monitor system memory usages 
  // dynamically.
  // Note: Before calling the following four functions, the user need to call 
  // PetscMemorySetGetMaximumUsage() immediately after PetscInitialize.
  // =================================================================
  // -----------------------------------------------------------------
  // 1. Print the Maximum memory used for the program. The usage is
  //    reduced to CPU 0 by MPI_SUM.
  // -----------------------------------------------------------------
  void print_MaxMemUsage();

  // ----------------------------------------------------------------
  // 2. Print the Current memory used for the program. The usage is 
  //    reduced to CPU 0 by MPI_SUM.  
  // ----------------------------------------------------------------
  void print_CurMemUsage();

  // ----------------------------------------------------------------
  // 3. Print the Maximum space PETSc has allocated. This function 
  //    should be used with the command line argument -malloc
  // ----------------------------------------------------------------
  void print_MaxMallocUsage();

  // ----------------------------------------------------------------
  // 4. Print the Current space PETSc has allocated. This function 
  //    should be used with the command line argument -malloc
  // ----------------------------------------------------------------
  void print_CurMallocUsage();

  // ================================================================
  // The following are system functions that access the system info.
  // ================================================================
  // 1. get_time: return the present time as HH:MM:SS
  // ----------------------------------------------------------------
  std::string get_time();

  // ----------------------------------------------------------------
  // 2. get_date: return the present date as YYYY/MM/DD
  // ----------------------------------------------------------------
  std::string get_date();

  // ----------------------------------------------------------------
  // 3. Structure that holds information about memory usage in kB.
  //    Used by get_memory_stats() function. See man 5 proc entry /
  //    status for details
  // ----------------------------------------------------------------
  struct MemoryStats
  {
    unsigned long int VmPeak; // peak virtual memory size in kB
    unsigned long int VmSize; // current virtual memory size in kB
    unsigned long int VmHWM; // peak resident memory size in kB
    unsigned long int VmRSS; // current resident memory size in kB
  };

  // ----------------------------------------------------------------
  // 4. get_memory_stats(): fills the MemoryStats structure with info
  //    about the memory consumption of this process. Only for Linux.
  // ----------------------------------------------------------------
  void get_memory_stats( MemoryStats &stats );

  // ================================================================
  // The folowing are options-get functions that read command-line
  // argument. They are really wrappers of the PetscOptionsGetXXXXX
  // functions.
  // ================================================================
  inline void GetOptionReal( const char * const &name, double &outdata )
  {
#if PETSC_VERSION_LT(3,7,0)
    PetscOptionsGetReal(PETSC_NULL, name, &outdata, PETSC_NULL);
#else
    PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, name, &outdata, PETSC_NULL);
#endif
  }

  inline void GetOptionInt( const char * const &name, int &outdata )
  {
#if PETSC_VERSION_LT(3,7,0)
    PetscOptionsGetInt(PETSC_NULL, name, &outdata, PETSC_NULL);
#else
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, name, &outdata, PETSC_NULL);
#endif
  }

  inline void GetOptionBool( const char * const &name, bool &outdata )
  {
    PetscBool pdata = PETSC_FALSE;
    PetscBool flg;
#if PETSC_VERSION_LT(3,7,0)
    PetscOptionsGetBool(PETSC_NULL, name, &pdata, &flg);
#else
    PetscOptionsGetBool(PETSC_NULL, PETSC_NULL, name, &pdata, &flg);
#endif
    if(flg)
    {
      if(pdata == PETSC_FALSE)
        outdata = false;
      else
        outdata = true;
    }
  }

  inline void GetOptionString( const char * const &name, std::string &outdata )
  {
    PetscBool flg;
    char char_outdata[PETSC_MAX_PATH_LEN];
#if PETSC_VERSION_LT(3,7,0)
    PetscOptionsGetString(PETSC_NULL, name, char_outdata, PETSC_MAX_PATH_LEN, &flg);
#else
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, name, char_outdata, PETSC_MAX_PATH_LEN, &flg);
#endif
    if(flg) outdata = char_outdata;
  }

  // ----------------------------------------------------------------
  // Check if a file exists. If a file cannot be found, throw an error
  // message and exit code
  // ----------------------------------------------------------------
  inline bool file_exist( const std::string &fName )
  {
    if( FILE *ff = fopen(fName.c_str(), "r") )
    {
      fclose(ff);
      return true;
    }
    else return false;
  }

  inline void file_check( const std::string &fName )
  {
    print_fatal_if( !file_exist(fName), 
        "Error: The file %s does not exist. Job is killed. \n", fName.c_str());
  }

  // ================================================================
  // SYS_T::Timer class defines a timer tool that one can use to
  // measure the time spent on events
  // ================================================================
  class Timer
  {
    public:
      Timer();
      ~Timer();

      void Start();
      void Stop();
      void Reset();

      double get_sec() const;

    private:
      clock_t startedAt;
      clock_t stoppedAt;
  };

}

#endif
