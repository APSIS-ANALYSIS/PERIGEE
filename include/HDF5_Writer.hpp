#ifndef HDF5_WRITER_HPP
#define HDF5_WRITER_HPP
// ==================================================================
// HDF5_Writer.hpp
// 
// Description:
// This is a collection of functions for writing int/double scalar/
// vector/tensor into HDF5 files.
// 
// The users should generate a file_id first by calling H5Fcreate or
// H5Fopen. After the output job is done by calling the HDF5_Writer member
// functions, the user should call H5Fclose.
// 
// A typical usage is:
// 
// hid_t file_id = H5Fcreate(name_of_h5_file, mode, H5P_DEFAULT, H5P_DEFAULT)
// or 
// hid_t file_id = H5Fopen(name_of_h5_file, mode, H5P_DEFAULT)
//
// call HDF5_Writer::functions
// 
// H5Fclose(file_id)
// 
// There are two modes for the H5Fcreate:
// H5F_ACC_TRUNC: if the file exists, the content will be deleted and
//                write new content
// H5F_ACC_EXCL: if the file exists, the open will fail
// 
// and two modes for the H5Fopen
// H5F_ACC_RDONLY: the application will read only the file
// H5F_ACC_RDWR: the application will read and write the file.
//
// For more details, see 
//      https://support.hdfgroup.org/HDF5/Tutor/crtfile.html
// 
// Author: Ju Liu
// Date: Oct. 23 2013
// ==================================================================
#include <string>
#include <cstdlib>
#include <climits>
#include "Tensor2_3D.hpp"
#include "hdf5.h"

class HDF5_Writer
{
  public:
    // --------------------------------------------------------------
    // constructor:
    // read in file_id. The user should call H5Fcreate to generate a 
    // valid file_id. After the writing, the user is also responsible 
    // to call H5Fclose to close the h5 file.
    // --------------------------------------------------------------
    HDF5_Writer( const hid_t &in_file_id ) : file_id(in_file_id) {}
    
    virtual ~HDF5_Writer() = default;
  
    // --------------------------------------------------------------
    // Scalar writer 
    // --------------------------------------------------------------
    // --- Write an int (H5T_NATIVE_INT) scalar as /data_name
    void write_intScalar( const char * const &data_name, 
        const int &value ) const; 

    // --- Write an int (H5T_NATIVE_INT) scalar as /group_id/data_name
    void write_intScalar( const hid_t &group_id, 
        const char * const &data_name, const int &value ) const; 

    // --- Write a double (H5T_NATIVE_DOUBLE) scalar as /file_id/group_id
    void write_doubleScalar( const hid_t &group_id, 
        const char * const &data_name, const double &value ) const; 
    
    // --- Write a double (H5T_NATIVE_DOUBLE) scalar as /file_id
    void write_doubleScalar( const char * const &data_name, 
        const double &value ) const; 
    
    // --- Write an int64_t (H5T_NATIVE_LLONG) scalar as /data_name
    void write_int64Scalar( const char * const &data_name, 
        const int64_t &value ) const; 
   
    // --- Write an unsigned int (H5T_NATIVE_UINT) scalar at 
    //     /group_id/data_name
    void write_uintScalar( const hid_t &group_id,
        const char * const &data_name, const unsigned int &value ) const;

    // --------------------------------------------------------------
    // Array writer
    //      32-bit Integer writer 
    // --------------------------------------------------------------
    // --- write an int array with given length and saved as 
    //     /group/data_name
    void write_intVector( const hid_t &group_id,
        const char * const &data_name, const int * const &value, 
        const int &length ) const;
    
    // --- write an int array with given length and saved as /data_name
    void write_intVector( const char * const &data_name, 
        const int * const &value, const int &length ) const;

    // --------------------------------------------------------------
    // Array writer
    //      64-bit Integer writer 
    // --------------------------------------------------------------
    // --- write an int array with given length and saved as 
    //     /group/data_name
    void write_int64Vector( const hid_t &group_id, 
        const char * const &data_name, 
        const int64_t * const &value, const int64_t &length ) const;
    
    // --- write an int array with given length and saved as /data_name
    void write_int64Vector( const char * const &data_name, 
        const int64_t * const &value, 
        const int64_t &length ) const;

    // --------------------------------------------------------------
    // Array writer
    //      double array writer
    // --------------------------------------------------------------
    // --- write a double array with given length and save as 
    //     /group/data_name
    void write_doubleVector( const hid_t &group_id, 
        const char * const &data_name,
        const double * const &value, const int &length ) const;

    // --- write a double array with given length at /data_name
    void write_doubleVector( const char * const &data_name,
        const double * const &value, const int &length ) const;

    // --------------------------------------------------------------
    // std::vector
    // --------------------------------------------------------------
    // --- write an int vector at /group_id/data_name
    void write_intVector( const hid_t &group_id,
        const char * const &data_name, const std::vector<int> &value ) const;

    // --- write an int vector at /data_name
    void write_intVector( const char * const &data_name, 
        const std::vector<int> &value ) const;
    
    // --- write an unsigned int vector at /group_id/data_name
    void write_uintVector( const hid_t &group_id,
        const char * const &data_name, 
        const std::vector<unsigned int> &value ) const;

    // --- write a double vector at /group_id/data_name
    void write_doubleVector( const hid_t &group_id,
        const char * const &data_name, const std::vector<double> &value ) const;

    // --- write a double vector at /data_name
    void write_doubleVector( const char * const &data_name, 
        const std::vector<double> &value ) const;

    // --------------------------------------------------------------
    // Vector_3
    // --------------------------------------------------------------
    void write_Vector_3( const hid_t &group_id, const char * const &data_name,
        const Vector_3 &value ) const;

    void write_Vector_3( const char * const &data_name, const Vector_3 &value ) const;

    // --------------------------------------------------------------
    // Tensor2_3D
    // --------------------------------------------------------------
    void write_Tensor2_3D( const hid_t &group_id, const char * const &data_name,
        const Tensor2_3D &value ) const;

    void write_Tensor2_3D( const char * const &data_name, const Tensor2_3D &value ) const;

    // --------------------------------------------------------------
    // Matrix writer
    // --------------------------------------------------------------
    void write_intMatrix( const hid_t &group_id,
        const char * const &data_name, const std::vector<int> &value,
        const int &row_num, const int &col_num ) const;

    void write_doubleMatrix( const hid_t &group_id, 
        const char * const &data_name, const std::vector<double> &value,
        const int &row_num, const int &col_num ) const;

    // --------------------------------------------------------------
    // String
    // --------------------------------------------------------------
    void write_string( const char * const &data_name, 
        const std::string &string_input ) const;

    void write_string( const hid_t &group_id,
        const char * const &data_name, 
        const std::string &string_input ) const;

  private:
    const hid_t file_id;

    void check_error(const herr_t &status, const char * const &funname ) const
    {
      if(status < 0) 
      {
        std::cerr<<"Error: HDF5_Writer::"<<funname<<std::endl;
        exit( EXIT_FAILURE );
      }
    }

    void write_string_impl(hid_t location_id, const char * const &data_name, 
        const std::string& string_input ) const;
};

#endif
