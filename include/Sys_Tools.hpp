#ifndef SYS_TOOLS_HPP
#define SYS_TOOLS_HPP
// ============================================================================
// Sys_Tools.hpp
// ----------------------------------------------------------------------------
// The SYS_T namespace contains a suite of tools at the system level.
// ============================================================================
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <sys/stat.h>
#include "petsc.h"
#ifdef USE_OPENMP
#include "omp.h"
#endif
#ifdef _OPENMP
#define PERIGEE_OMP_PARALLEL_FOR _Pragma("omp parallel for")
#define PERIGEE_OMP_PARALLEL _Pragma("omp parallel")
#define PERIGEE_OMP_FOR _Pragma("omp for")
#define PERIGEE_OMP_CRITICAL _Pragma("omp critical")
#define PERIGEE_OMP_SINGLE _Pragma("omp single")
#else
#define PERIGEE_OMP_PARALLEL_FOR
#define PERIGEE_OMP_PARALLEL
#define PERIGEE_OMP_FOR
#define PERIGEE_OMP_CRITICAL
#define PERIGEE_OMP_SINGLE
#endif

#define PETSC_SILENCE_DEPRECATION_WARNINGS_3_19_0

  // ================================================================
  // The following are used for backward compatibility like PetscDefined(USE_DEBUG).
  // ================================================================
  // Versions >= 3.14.x : PetscDefined(USE_DEBUG) is used to determine whether it is debug mode;
  //           < 3.14.x : defined(PETSC_USE_DEBUG) is used to determine whether it is debug mode.
#if PETSC_VERSION_LT(3,14,6)
  #define PETSC_DEFINED(def) defined(PETSC_ ## def)
#else
  #define PETSC_DEFINED(def) PetscDefined(def)
#endif

  // ================================================================
  // The following are used for ASSERT.
  // ================================================================
  // In debug mode, ASSERT is called to determine a "cond" condition.
#if PETSC_DEFINED(USE_DEBUG)
  #define ASSERT(cond, message, ...) SYS_T::print_fatal_if_not(cond, message, ##__VA_ARGS__)
#else
  #define ASSERT(cond, ...) ((void)0)
#endif

namespace SYS_T
{
  // Return the OS environmental variable
  inline std::string get_Env_Var( const std::string &key )
  {
    const char * val = std::getenv( key.c_str() );
    return val == nullptr ? std::string("") : std::string(val);
  }

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

  // ================================================================
  // The following are used in processor.
  // ================================================================
  // 1. Synchronized print on screen. Here we use petsc functions:
  //    PetscSynchronizedPrintf() and PetscSynchronizedFlush().
  //    The output should be lised from proc 0, 1, ... to n in sequence.
  inline void synPrint(const std::string &output, const int &cpu_rank)
  {
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Proc %d: ", cpu_rank);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", output.c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  }
  inline void synPrint(const char * const &output, const int &cpu_rank)
  {
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Proc %d: ", cpu_rank);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", output);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  }
  inline void synPrint(const std::string &output)
  {
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", output.c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  }
  inline void synPrint(const char * const &output)
  {
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", output);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  }

  // 2. print from processor 0, other preprocessors are ignored.
  inline void commPrint(const char output[], ...)
  {
    int mpi_flag {-1};
    MPI_Initialized(&mpi_flag);
    if( mpi_flag )
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
    else
    {
#ifdef _OPENMP
      if( !omp_get_thread_num() )
      {
        va_list Argp;
        va_start(Argp, output);
        vfprintf(stdout, output, Argp);
        va_end(Argp);
      }
      #pragma omp barrier
#else
      va_list Argp;
      va_start(Argp, output);
      vfprintf(stdout, output, Argp);
      va_end(Argp);
#endif
    }
  }

  // 3. Print the data name and its value on screen using PetscPrintf
  //    and PETSC_COMM_WORLD. This is particularly designed to print
  //    command line arguments.
  inline void cmdPrint(const char * const &dataname, const int &datavalue)
  {std::ostringstream ss; ss<<dataname<<" "<<datavalue<<"\n"; PetscPrintf(PETSC_COMM_WORLD, "%s", ss.str().c_str());}
  inline void cmdPrint(const char * const &dataname, const double &datavalue)
  {std::ostringstream ss; ss<<dataname<<" "<<datavalue<<"\n"; PetscPrintf(PETSC_COMM_WORLD, "%s", ss.str().c_str());}
  inline void cmdPrint(const char * const &dataname, const std::string &datavalue)
  {std::ostringstream ss; ss<<dataname<<" "<<datavalue<<"\n"; PetscPrintf(PETSC_COMM_WORLD, "%s", ss.str().c_str());}

  // 4. Print fatal error message and terminate the MPI process
  inline void print_fatal( const char output[], ... )
  { 
    int mpi_flag {-1};
    MPI_Initialized(&mpi_flag);
    if (mpi_flag)
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
    else
    {
#ifdef _OPENMP
      if( !omp_get_thread_num() )
      {
        va_list Argp;
        va_start(Argp, output);
        vfprintf (stderr, output, Argp);
        va_end(Argp);

        exit( EXIT_FAILURE );
      }
      else exit( EXIT_FAILURE );
#else
      va_list Argp;
      va_start(Argp, output);
      vfprintf (stderr, output, Argp);
      va_end(Argp);

      exit( EXIT_FAILURE );
#endif
    }
  }

  inline void print_fatal_if( bool a, const char output[], ... )
  {
    if( a )
    {
      int mpi_flag {-1};
      MPI_Initialized(&mpi_flag);
      if (mpi_flag)
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
      else
      {
#ifdef _OPENMP
        if( !omp_get_thread_num() )
        {
          va_list Argp;
          va_start(Argp, output);
          vfprintf (stderr, output, Argp);
          va_end(Argp);

          exit( EXIT_FAILURE );
        }
        else exit( EXIT_FAILURE );
#else
        va_list Argp;
        va_start(Argp, output);
        vfprintf (stderr, output, Argp);
        va_end(Argp);

        exit( EXIT_FAILURE );
#endif
      }      
    }
  }

  inline void print_fatal_if_not( bool a, const char output[], ... )
  {
    if( !a )
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

  // 4. print the number of threads used in openmp
  inline void print_omp_info()
  {
#ifdef _OPENMP
    PERIGEE_OMP_PARALLEL
    {
      PERIGEE_OMP_SINGLE
      {
        std::cout<<"The number of threads used: "<<omp_get_num_threads()<<", and ";
        std::cout<<"the number of processors on the machine: ";
        std::cout<<omp_get_num_procs()<<".\n";
      }
    }
#else
    std::cout<<"OpenMP is not invoked.\n";
#endif
  }

  // 5. set the number of threads used in openmp
  inline void set_omp_num_threads()
  {
#ifdef _OPENMP
    omp_set_num_threads( omp_get_num_procs() );
#endif
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

  inline void InsertFileYAML( const std::string &filename, const bool &require )
  {
#if PETSC_VERSION_GE(3,15,0)
    if( require )
    {
      file_check( filename );
      commPrint("Status: loading YAML file %s from the disk.\n", filename.c_str());
      PetscOptionsInsertFileYAML(PETSC_COMM_WORLD, NULL, filename.c_str(), PETSC_TRUE);
    }
    else
    {
      if(!file_exist(filename))
        commPrint("Warning: the YAML file %s does not exist, and the command line arguments are loaded.\n", filename.c_str());

      PetscOptionsInsertFileYAML(PETSC_COMM_WORLD, NULL, filename.c_str(), PETSC_FALSE);
    }
#else
    commPrint("Warning: YAML is unsupported in this PETSc.\n");
#endif
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

      void Reset() { startedAt = 0; stoppedAt = 0; }

#ifdef _OPENMP
      void Start() { startedAt = omp_get_wtime(); }
      void Stop()  { stoppedAt = omp_get_wtime(); }
      double get_sec() const
      {
        return (stoppedAt - startedAt);
      }
#else
      void Start() { startedAt = clock(); }
      void Stop()  { stoppedAt = clock(); }
      double get_sec() const
      {
        return (double)(stoppedAt - startedAt)/(double)CLOCKS_PER_SEC;
      }
#endif

    private:
#ifdef _OPENMP
      double startedAt, stoppedAt;
#else
      clock_t startedAt, stoppedAt;
#endif
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

  inline void print_system_info()
  {
    commPrint("Date: %s and time: %s.\n", get_date().c_str(), get_time().c_str());
    commPrint("Machine: %s \n", get_Env_Var("MACHINE_NAME").c_str());
    commPrint("User: %s \n", get_Env_Var("USER").c_str());
    commPrint("The sizes of int, double, and long double are %zu byte, %zu byte, and %zu byte, resp.\n", sizeof(int), sizeof(double), sizeof(long double) );
  }

}

#endif
