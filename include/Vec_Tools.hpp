#ifndef VEC_TOOLS_HPP
#define VEC_TOOLS_HPP
// ============================================================================
// Vec_Tools.hpp
// ----------------------------------------------------------------------------
// VEC_T namespace contains a suite of functions for the manipulation of the 
// std::vector object.
// ============================================================================
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

namespace VEC_T
{
  // --------------------------------------------------------------------------
  // ! print
  //   print int or double vector on screen. By default, the separation of the 
  //   vector entries is \t. User may switch it to \n.
  // --------------------------------------------------------------------------
  template<typename T> void print( const std::vector<T> &vec, const char &sep = '\t' )
  {
    for( auto it = vec.begin(); it != vec.end(); ++it )
      std::cout<<std::setprecision(16)<<*it<<sep;
    std::cout<<'\n';
  }
 
  // --------------------------------------------------------------------------
  // ! print
  //   print the vector with given precision pres 
  // --------------------------------------------------------------------------
  template<typename T> void print( const std::vector<T> &vec, 
      const unsigned int &pres, const char &sep = '\t' )
  {
    const std::streamsize ss = std::cout.precision();
    for( auto it = vec.begin(); it != vec.end(); ++it )
      std::cout<<std::setprecision(pres)<<*it<<sep;
    std::cout<<'\n';
    std::cout.precision(ss);
  }
  
  // --------------------------------------------------------------------------
  // ! print
  //   print the vector to a file with given file name
  // --------------------------------------------------------------------------
  template<typename T> void print( const std::vector<T> &vec,
      const std::string &file_name, const char &sep = '\t'  )
  {
    std::ofstream efile(file_name.c_str(), std::ofstream::out | std::ofstream::trunc );
    for(auto it = vec.begin(); it!= vec.end(); ++it) efile<<*it<<sep;
    
    efile.close();
  }

  // --------------------------------------------------------------------------
  // ! get_size
  //   return the length of a std::vector in signed int type (rather than 
  //   unsigned int).
  // --------------------------------------------------------------------------
  template<typename T> int get_size( const std::vector<T> &vec )
  {
    return static_cast<int>(vec.size());
  }

  // --------------------------------------------------------------------------
  // ! is_equal
  // determine if two vector object are identical up to a tolerance.
  // --------------------------------------------------------------------------
  template<typename T> bool is_equal( const std::vector<T> &a, 
      const std::vector<T> &b, double tol = 1.0e-12 )
  {
    if( a.size() != b.size() ) return false;
    for(unsigned int ii=0; ii<a.size(); ++ii)
    {
      if( std::abs(a[ii]-b[ii]) >= tol ) return false;
    }
    return true;
  }
  
  // --------------------------------------------------------------------------
  // ! fillArray
  //   fill the array data into a vector, with array length len
  // --------------------------------------------------------------------------
  template<typename T> std::vector<T> fillArray( const T * const &input, int len )
  {
    std::vector<T> vec(len);
    for(int ii=0; ii<len; ++ii) vec[ii] = input[ii];
    
    vec.shrink_to_fit();
    return vec;
  }

  // --------------------------------------------------------------------------
  // ! max
  //   return the maximum value in the vector
  // --------------------------------------------------------------------------
  template<typename T> T max( const std::vector<T> &vec )
  {
    return *std::max_element( vec.begin(), vec.end() );
  }

  // --------------------------------------------------------------------------
  // ! min
  //   return the minimum value in the vector
  // --------------------------------------------------------------------------
  template<typename T> T min( const std::vector<T> &vec )
  {
    return *std::min_element( vec.begin(), vec.end() );
  }

  // --------------------------------------------------------------------------
  // ! sum
  //   return the summation of all entries in the vector
  // --------------------------------------------------------------------------
  template<typename T> T sum( std::vector<T> &vec )
  {
    T counter = 0;
    for( const auto &entry : vec ) counter += entry;
    return counter;
  }

  // --------------------------------------------------------------------------
  // ! insert_end
  //   insert vec_b at the end of vec_a
  // --------------------------------------------------------------------------
  template<typename T> void insert_end( std::vector<T> &vec_a,
      const std::vector<T> &vec_b )
  {
    vec_a.insert(vec_a.end(), vec_b.begin(), vec_b.end());
  } 

  // --------------------------------------------------------------------------
  // ! clean
  //   clean the allocation of a vector
  // --------------------------------------------------------------------------
  template<typename T> void clean( std::vector<T> &vec ) 
  {std::vector<T>().swap(vec);} 

  // --------------------------------------------------------------------------
  // ! sort_unique_resize
  //   \para vec : the vector to be sorted and resized.
  //   sort the vector first, 
  //   then delete repeated items and resize the vector
  //   e.g. vec = [0,3,4,1,3,4,5,7,8,9,3]
  //   will become
  //        vec = [0,1,3,4,5,7,8,9].
  // --------------------------------------------------------------------------
  template<typename T> void sort_unique_resize( std::vector<T> &vec )
  {
    sort(vec.begin(), vec.end());
    vec.resize( unique(vec.begin(), vec.end()) - vec.begin() );
  }

  // --------------------------------------------------------------------------
  // ! set_diff
  //   input: vec_a and vec_b
  //   output : the set difference of vec_a - vec_b, that is the value belonging
  //   to vec_a, but not vec_b.
  //   e.g. vec_a = [ 5, 10, 15, 20, 10, 25 ]; vec_b = [ 10, 20, 30, 50, 30, 40 ];
  //   output is [ 5, 15, 25 ].
  // --------------------------------------------------------------------------
  template<typename T> std::vector<T> set_diff( const std::vector<T> &vec_a,
      const std::vector<T> &vec_b )
  {
    auto temp_a = vec_a, temp_b = vec_b;

    sort_unique_resize( temp_a );
    sort_unique_resize( temp_b );

    auto output = temp_a;
    auto it = std::set_difference( temp_a.begin(), temp_a.end(), temp_b.begin(), temp_b.end(), output.begin() );
    output.resize( it - output.begin() );
    return output;
  }

  // --------------------------------------------------------------------------
  // ! is_invec
  //   determine if a given value val is in the vector vec (return true),
  //   or not (return false).
  // --------------------------------------------------------------------------
  template<typename T> bool is_invec( const std::vector<T> &vec, const T &val )
  {
    return find(vec.begin(), vec.end(), val) != vec.end();
  }

  // --------------------------------------------------------------------------
  // ! intersection
  //   return the intersection set of two given vectors
  // --------------------------------------------------------------------------
  template<typename T> std::vector<T> intersection( const std::vector<T> &vec_a,
     const std::vector<T> &vec_b )
  {
    std::vector<T> out {};
    for( const auto &val : vec_a )
      if( VEC_T::is_invec(vec_b, val) ) out.push_back(val);

    return out;
  }

  // --------------------------------------------------------------------------
  // ! get_pos
  //   find the position (pos) of the given val such that 
  //                 vec[pos] = val
  //   if there are multiple same values in vector, it will return the first 
  //   appearance of the given val.
  //   if the val is not found, return -1.
  // --------------------------------------------------------------------------
  template<typename T> int get_pos( const std::vector<T> &vec, const T &val )
  {
    const auto it = find(vec.begin(), vec.end(), val);
    if( it == vec.end() ) return -1;
    else return it - vec.begin();
  }

  // --------------------------------------------------------------------------
  // ! cast_to_unsigned_int
  //   Convert a std::vector<T> to std::vector<unsigned int>.
  // --------------------------------------------------------------------------
  template<typename T> std::vector<unsigned int> cast_to_unsigned_int( const std::vector<T> &vec )
  {
    std::vector<unsigned int> output( vec.size() );

    for(unsigned int ii=0; ii<vec.size(); ++ii)
      output[ii] = static_cast<unsigned int>( vec[ii] );
    
    return output;
  }

  // --------------------------------------------------------------------------
  // ! write_txt
  //   write a vector to disk in .txt format.
  //   \para file_name : file_name.vec.txt will be the name for the file
  //   \para vec : the vector to be written
  //   Note: this file can be read by matlab and do various operations in
  //   matlab.
  // --------------------------------------------------------------------------
  template<typename T> void write_txt( const std::string &file_name,
      const std::vector<T> &vec, const bool &wIdx=true, const unsigned int &pres=6 )
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

  // --------------------------------------------------------------------------
  // ! erase_pos
  //   erase the element stored in the vector vec at the position at begin() +
  //   pos
  // --------------------------------------------------------------------------
  template<typename T> void erase_pos( std::vector<T> &vec, int pos )
  {
    vec.erase( vec.begin() + pos );
  }

}

#endif
