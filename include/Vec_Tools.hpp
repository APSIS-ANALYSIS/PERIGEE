#ifndef VEC_TOOLS_HPP
#define VEC_TOOLS_HPP
// ==================================================================
// Vec_Tools.hpp
// ------------------------------------------------------------------
// VEC_T namespace contains a suite of frequently used function for
// the std::vector object.
// ==================================================================
#include <iomanip>
#include <algorithm>
#include "Sys_Tools.hpp"
#include "HDF5_Writer.hpp"
#include "HDF5_Reader.hpp"

namespace VEC_T
{
  // print int / double vector on screen
  template<typename T> void print( const std::vector<T> &vec )
  {
    for( auto it = vec.begin(); it != vec.end(); ++it )
      std::cout<<std::setprecision(16)<<*it<<'\t';
    std::cout<<'\n';
  }
 
  
  // Print the vector with given precision pres 
  template<typename T> void print( const std::vector<T> &vec, 
      const unsigned int pres )
  {
    const std::streamsize ss = std::cout.precision();
    for( auto it = vec.begin(); it != vec.end(); ++it )
      std::cout<<std::setprecision(pres)<<*it<<'\t';
    std::cout<<'\n';
    std::cout.precision(ss);
  }
  
  
  // Print the vector to a file with given file name
  template<typename T> void print( const std::vector<T> &vec,
      const std::string &file_name )
  {
    std::ofstream efile(file_name.c_str(), std::ofstream::out | std::ofstream::trunc );
    for(auto it = vec.begin(); it!= vec.end(); ++it) efile<<*it<<'\t';
    
    efile.close();
  }


  // trim the capacity of vector
  template<typename T> void shrink2fit( std::vector<T> &vec )
  {
    std::vector<T>(vec.begin(), vec.end()).swap(vec);
  }


  // fill array data into a vector
  template<typename T> void fillArray( std::vector<T> &vec,
      const T * const &input, const int &len )
  {
    vec.clear(); vec.resize(len);
    for(int ii=0; ii<len; ++ii) vec[ii] = input[ii];
    
    shrink2fit(vec);
  }


  // insert vec_b at the end of vec_a
  template<typename T> void insert_end( std::vector<T> &vec_a,
      const std::vector<T> &vec_b )
  {
    vec_a.insert(vec_a.end(), vec_b.begin(), vec_b.end());
  } 

  
  // clean the allocation of a vector
  template<typename T> void clean( std::vector<T> &vec )
  {std::vector<T>().swap(vec);} 


  // generate a double vector with random entries
  // default length is 1, random double in [min=0, max=1]
  inline void gen_random_double( std::vector<double> &vec,
      const unsigned int &vec_len = 1, const int &min = 0, 
      const int &max = 1 )
  {
    vec.resize(vec_len);
    srand(time(NULL));
    for(unsigned int ii=0; ii<vec_len; ++ii) 
      vec[ii] = SYS_T::gen_randomD_closed(min, max);
  }

  
  // generate a integer vector with random entries,
  // default length is 1, random int in [min=0, max=1]
  inline void gen_random_int( std::vector<int> &vec,
      const unsigned int &vec_len = 1, const int &min = 0,
      const int &max = 1 )
  {
    vec.resize(vec_len);
    srand(time(NULL));
    for(unsigned int ii=0; ii<vec_len; ++ii) 
      vec[ii] = SYS_T::gen_randomI_closed(min, max);
  }


  // inner product of two vectors with given length len
  inline double vDot( const std::vector<double> &vec1, 
      const std::vector<double> &vec2, const unsigned int &len)
  {
    double result = 0;
    for(unsigned int ii=0; ii<len; ++ii)
      result += vec1[ii] * vec2[ii];
    return result;
  }


  // ----------------------------------------------------------------
  // ! sort_unique_resize
  //   \para vec : the vector to be sorted and resized.
  //   sort the vector first, 
  //   then delete repeated items and resize the vector
  //   e.g. vec = [0,3,4,1,3,4,5,7,8,9,3]
  //   will become
  //        vec = [0,1,3,4,5,7,8,9].
  // ----------------------------------------------------------------
  template<typename T> void sort_unique_resize( std::vector<T> &vec )
  {
    sort(vec.begin(), vec.end());
    auto ite = unique(vec.begin(), vec.end());
    vec.resize( ite - vec.begin() );
  }

  
  // ---------------------------------------------------------------- 
  // ! is_invec
  //   determine if a given value val is in the vector vec.
  // ---------------------------------------------------------------- 
  template<typename T> bool is_invec( const std::vector<T> &vec, 
      const T &val)
  {
    auto it = find(vec.begin(), vec.end(), val);
    return it != vec.end();
  }


  // ----------------------------------------------------------------
  // ! is_equal
  //   determine if two given vectors are the same
  // ----------------------------------------------------------------
  template<typename T> bool is_equal( const std::vector<T> &vec_a,
      const std::vector<T> &vec_b )
  {
    if(vec_a.size() != vec_b.size() ) return false;
    for(unsigned int ii=0; ii<vec_a.size(); ++ii)
    {
      if(vec_a[ii] != vec_b[ii]) return false;
    }
    return true;
  }


  // ----------------------------------------------------------------
  // ! get_pos
  //   find the position (pos) of the given val such that 
  //                 vec[pos] = val
  //   if there are multiple same values in vector, it will return the
  //   first appearance of the givne val.
  //   if the val is not found, return -1.
  // ----------------------------------------------------------------
  template<typename T> int get_pos( const std::vector<T> &vec,
      const T &val )
  {
    const auto it = find(vec.begin(), vec.end(), val);
    if( it == vec.end() ) return -1;
    else return it - vec.begin();
  }


  // -----------------------------------------------------------------
  // ! write_txt
  //   write a vector to disk in .txt format.
  //   \para file_name : file_name.vec.txt will be the name for the file
  //   \para vec : the vector to be written
  //   Note: this file can be read by matlab and do various operations in
  //   matlab.
  // -----------------------------------------------------------------
  template<typename T> void write_txt( const char * const &file_name,
      const std::vector<T> &vec, const bool &wIdx=true,
      const unsigned int &pres=6 )
  {
    const std::streamsize ss = std::cout.precision();
    std::string fname(file_name);
    fname.append(".vec.txt");
    std::ofstream vfile;
    vfile.open( fname.c_str(), std::ofstream::out | std::ofstream::trunc );
    if(wIdx)
    {
      for(unsigned int ii=0; ii<vec.size(); ++ii)
        vfile <<std::setprecision(pres) << ii <<'\t'<< vec[ii] << std::endl;
    }
    else
    {
      for(unsigned int ii=0; ii<vec.size(); ++ii)
        vfile <<std::setprecision(pres) << vec[ii] << std::endl;
    }
    vfile.close();
    std::cout.precision(ss);
  }


  // -----------------------------------------------------------------
  // ! write_int_h5
  //   Write an int vector to disk as HDF5 file.
  //   \para file_name : the input that will name the .h5 file
  //   \para dataname : the name of the dataset in the file
  //   \para value : the vector that is to be written.
  // -----------------------------------------------------------------
  void write_int_h5( const char * const &file_name, 
      const char * const &dataname,
      const std::vector<int> &value );


  // -----------------------------------------------------------------
  // ! read_int_h5
  //   Read an int vector from a HDF5 file on disk
  //   \para file_name : the file's name
  //   \para groupname : the group name of the dataset
  //   \para dataname : the dataset's name in the file
  //   \para value : the vector that will be filled in.
  // -----------------------------------------------------------------
  void read_int_h5( const char * const &file_name, 
      const char * const &groupname,
      const char * const &dataname,
      std::vector<int> &value );
}

#endif
