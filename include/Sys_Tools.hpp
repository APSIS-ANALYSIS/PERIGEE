#ifndef SYS_TOOLS_HPP
#define SYS_TOOLS_HPP
// ============================================================================
// Sys_Tools.hpp
// ----------------------------------------------------------------------------
// The SYS_T namespace contains a suite of tools at the system level.
//
// Author: Ju Liu
// ============================================================================
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <iomanip>
#include <sys/stat.h>
#include "petsc.h"

namespace SYS_T
{
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

  // --------------------------------------------------------------------------
  // ! gen_partfile_name( baseName, rank )
  //   Generate a partition file's name (hdf5 file) in the form baseName_pxxxxx.h5.
  // --------------------------------------------------------------------------
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

  // --------------------------------------------------------------------------
  // ! gen_capfile_name
  //   Generate a file (usually for cap surfaces) in the form baseName_xxx.vtp(vtu)
  // --------------------------------------------------------------------------
  inline std::string gen_capfile_name( const std::string &baseName,
      const int &index, const std::string &filename )
  {
    std::ostringstream ss;
    ss<<baseName;

    if( index/10 == 0 ) ss<<"00";
    else if( index/100 == 0 ) ss<<"0";

    ss<<index<<filename;

    return ss.str();
  }

  // --------------------------------------------------------------------------
  // print size based on bytes
  // --------------------------------------------------------------------------
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

  template<typename T> void print_mem_size( const T &byte_size )
  {
    std::cout<<get_string_mem_size( byte_size );
  }
  
  inline std::string to_string( const int &aa )
  {
    std::ostringstream ss;
    ss<<aa;
    return ss.str();
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
    if( !get_MPI_rank() )
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

  // 4. Print fatal error message and terminate the MPI process
  inline void print_fatal( const char output[], ... )
  {
    if( !get_MPI_rank() )
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
      if( !get_MPI_rank() )
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

  // 5. Print message (without termination the code) under conditions
  inline void print_message_if( bool a, const char output[], ... )
  {
    if( a )
    {
      if( !get_MPI_rank() )
      {
        va_list Argp;
        va_start(Argp, output);
        (*PetscVFPrintf)(PETSC_STDOUT,output,Argp);
        va_end(Argp);
      }

      MPI_Barrier(PETSC_COMM_WORLD);
    }
  }

  // 6. Print exit message printers are used in terminating serial 
  //    program when the communicator for MPI is not available.
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

  // ================================================================
  // The following are system functions that access the system info.
  // ================================================================
  // 1. get_time: return the present time as HH:MM:SS
  // ----------------------------------------------------------------
  inline std::string get_time()
  {
    std::time_t  time1= std::time (0);
    std::tm     *time = std::localtime(&time1);

    std::ostringstream o;
    o << time->tm_hour << ":"
      << (time->tm_min < 10 ? "0" : "") << time->tm_min << ":"
      << (time->tm_sec < 10 ? "0" : "") << time->tm_sec;

    return o.str();
  }

  // ----------------------------------------------------------------
  // 2. get_date: return the present date as YYYY/MM/DD
  // ----------------------------------------------------------------
  inline std::string get_date()
  {
    std::time_t  time1= std::time (0);
    std::tm     *time = std::localtime(&time1);

    std::ostringstream o;
    o << time->tm_year + 1900 << "/"
      << time->tm_mon + 1 << "/"
      << time->tm_mday;

    return o.str();
  }

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
  inline void get_memory_stats( MemoryStats &stats )
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

  inline bool directory_exist( const std::string &dName )
  {
    if (dName.empty() || dName == "" || dName == "/0")
      return true;
    
    struct stat info;
    if (stat(dName.c_str(), &info) == 0)
      return true;

    return false;
  }

  inline void file_check( const std::string &fName )
  {
    print_fatal_if( !file_exist(fName), 
        "Error: The file %s does not exist. Job is killed. \n", fName.c_str());
  }

  // --------------------------------------------------------------------------
  // Execute a system call
  // --------------------------------------------------------------------------
  inline void execute( const char * const &command )
  {
    int sysret = system( command );
    print_fatal_if(sysret != 0, "Error: system call %s failed. \n", command);
  }

  // ================================================================
  // SYS_T::Timer class defines a timer tool that one can use to
  // measure the time spent on events
  // ================================================================
  class Timer
  {
    public:
      Timer() { startedAt = 0; stoppedAt = 0; }

      ~Timer() {};

      void Start() {startedAt = clock();}

      void Stop() {stoppedAt = clock();}

      void Reset() { startedAt = 0; stoppedAt = 0; }

      double get_sec() const
      {
        return (double)(stoppedAt - startedAt)/(double)CLOCKS_PER_SEC;
      }

    private:
      clock_t startedAt, stoppedAt;
  };

  // Print ASCII art text for the code
  inline void print_perigee_art()
  {
    commPrint("$$$$$$$\\  $$$$$$$$\\ $$$$$$$\\  $$$$$$\\  $$$$$$\\  $$$$$$$$\\ $$$$$$$$\\ \n");
    commPrint("$$  __$$\\ $$  _____|$$  __$$\\ \\_$$  _|$$  __$$\\ $$  _____|$$  _____| \n");
    commPrint("$$ |  $$ |$$ |      $$ |  $$ |  $$ |  $$ /  \\__|$$ |      $$ | \n");
    commPrint("$$$$$$$  |$$$$$\\    $$$$$$$  |  $$ |  $$ |$$$$\\ $$$$$\\    $$$$$\\ \n");
    commPrint("$$  ____/ $$  __|   $$  __$$<   $$ |  $$ |\\_$$ |$$  __|   $$  __| \n");
    commPrint("$$ |      $$ |      $$ |  $$ |  $$ |  $$ |  $$ |$$ |      $$ | \n");
    commPrint("$$ |      $$$$$$$$\\ $$ |  $$ |$$$$$$\\ \\$$$$$$  |$$$$$$$$\\ $$$$$$$$\\ \n");
    commPrint("\\__|      \\________|\\__|  \\__|\\______| \\______/ \\________|\\________| \n \n");
  }

  inline void print_sep_line()
  {
    commPrint("----------------------------------------------------------------------\n");
  }
  
  inline void print_sep_double_line()
  {
    commPrint("======================================================================\n");
  }
}

#endif
