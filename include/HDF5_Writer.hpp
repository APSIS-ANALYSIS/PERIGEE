#ifndef HDF5_WRITER_HPP
#define HDF5_WRITER_HPP
// ==================================================================
// HDF5_Writer.hpp
// 
// Description:
// This is a collection of functions for writing int/double scalar/
// vector/tensor into HDF5 files.
// 
// The users should provide the h5 file name and open/create mode when
// constructing HDF5_Writer. The HDF5 file will be closed in the destructor.
// 
// A typical usage is:
// 
// HDF5_Writer h5w(name_of_h5_file, mode)
//
// call h5w.functions
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
#include <cstdlib>
#include <climits>
#include <vector>
#include <array>
#include "Sys_Tools.hpp"
#include "HDF5_Group.hpp"
#include "hdf5.h"

class HDF5_Writer
{
  public:
    // --------------------------------------------------------------
    // constructor:
    // open/create a h5 file from file name and mode.
    // for create modes (H5F_ACC_TRUNC/H5F_ACC_EXCL), call H5Fcreate;
    // for open modes (H5F_ACC_RDONLY/H5F_ACC_RDWR), call H5Fopen.
    // default mode is H5F_ACC_TRUNC.
    // --------------------------------------------------------------
    HDF5_Writer( const std::string &in_file_name,
        const unsigned in_mode = H5F_ACC_TRUNC );
    
    virtual ~HDF5_Writer();

    hid_t get_file_id() const { return file_id; }
  
    // --------------------------------------------------------------
    // Scalar writer 
    // --------------------------------------------------------------
    // --- Write an int (H5T_NATIVE_INT) scalar as /data_name
    void write_intScalar( const char * data_name, int value ) const; 

    // --- Write an int (H5T_NATIVE_INT) scalar as /group_id/data_name
    void write_intScalar( hid_t group_id, const char * data_name, int value ) const; 
    void write_intScalar( const HDF5_Group &group, const char * data_name, int value ) const
    { write_intScalar( group.id(), data_name, value ); }

    // --- Write a double (H5T_NATIVE_DOUBLE) scalar as /file_id/group_id
    void write_doubleScalar( hid_t group_id, const char * data_name, 
        double value ) const; 
    void write_doubleScalar( const HDF5_Group &group, const char * data_name,
        double value ) const
    { write_doubleScalar( group.id(), data_name, value ); }
    
    // --- Write a double (H5T_NATIVE_DOUBLE) scalar as /file_id
    void write_doubleScalar( const char * data_name, double value ) const;
    
    // --- Write an int64_t (H5T_NATIVE_LLONG) scalar as /data_name
    void write_int64Scalar( const char * data_name, int64_t value ) const;
   
    // --- Write an unsigned int (H5T_NATIVE_UINT) scalar at 
    //     /group_id/data_name
    void write_uintScalar( hid_t group_id, const char * data_name, 
        unsigned int value ) const;
    void write_uintScalar( const HDF5_Group &group, const char * data_name,
        unsigned int value ) const
    { write_uintScalar( group.id(), data_name, value ); }

    // --------------------------------------------------------------
    // Array writer
    //      32-bit Integer writer 
    // --------------------------------------------------------------
    // --- write an int array with given length and saved as 
    //     /group/data_name
    void write_intVector( hid_t group_id, const char * data_name,
        const int * value, int length ) const;
    void write_intVector( const HDF5_Group &group, const char * data_name,
        const int * value, int length ) const
    { write_intVector( group.id(), data_name, value, length ); }
    
    // --- write an int array with given length and saved as /data_name
    void write_intVector( const char * data_name, const int * value,
        int length ) const;

    // --------------------------------------------------------------
    // Array writer
    //      64-bit Integer writer 
    // --------------------------------------------------------------
    // --- write an int array with given length and saved as 
    //     /group/data_name
    void write_int64Vector( hid_t group_id, const char * data_name, 
        const int64_t * value, int64_t length ) const;
    void write_int64Vector( const HDF5_Group &group, const char * data_name,
        const int64_t * value, int64_t length ) const
    { write_int64Vector( group.id(), data_name, value, length ); }
    
    // --- write an int array with given length and saved as /data_name
    void write_int64Vector( const char * data_name, 
        const int64_t * value, int64_t length ) const;

    // --------------------------------------------------------------
    // Array writer
    //      double array writer
    // --------------------------------------------------------------
    // --- write a double array with given length and save as 
    //     /group/data_name
    void write_doubleVector( hid_t group_id, 
        const char * data_name,
        const double * value, int length ) const;
    void write_doubleVector( const HDF5_Group &group, const char * data_name,
        const double * value, int length ) const
    { write_doubleVector( group.id(), data_name, value, length ); }

    // --- write a double array with given length at /data_name
    void write_doubleVector( const char * data_name,
        const double * value, int length ) const;

    // --------------------------------------------------------------
    // std::vector
    // --------------------------------------------------------------
    // --- write an int vector at /group_id/data_name
    void write_intVector( hid_t group_id,
        const char * data_name, const std::vector<int> &value ) const;
    void write_intVector( const HDF5_Group &group, const char * data_name,
        const std::vector<int> &value ) const
    { write_intVector( group.id(), data_name, value ); }

    // --- write an int vector at /data_name
    void write_intVector( const char * data_name, 
        const std::vector<int> &value ) const;
    
    // --- write an unsigned int vector at /group_id/data_name
    void write_uintVector( hid_t group_id, const char * data_name, 
        const std::vector<unsigned int> &value ) const;
    void write_uintVector( const HDF5_Group &group, const char * data_name,
        const std::vector<unsigned int> &value ) const
    { write_uintVector( group.id(), data_name, value ); }

    // --- write a double vector at /group_id/data_name
    void write_doubleVector( hid_t group_id, const char * data_name,
        const std::vector<double> &value ) const;
    void write_doubleVector( const HDF5_Group &group, const char * data_name,
        const std::vector<double> &value ) const
    { write_doubleVector( group.id(), data_name, value ); }

    // --- write a double vector at /data_name
    void write_doubleVector( const char * data_name, 
        const std::vector<double> &value ) const;

    // --------------------------------------------------------------
    // Vector_3
    // --------------------------------------------------------------
    void write_Vector_3( hid_t group_id, const char * data_name,
        const std::array<double, 3> &value ) const;
    void write_Vector_3( const HDF5_Group &group, const char * data_name,
        const std::array<double, 3> &value ) const
    { write_Vector_3( group.id(), data_name, value ); }

    void write_Vector_3( const char * data_name, 
        const std::array<double, 3> &value ) const;

    // --------------------------------------------------------------
    // Tensor2_3D
    // --------------------------------------------------------------
    void write_Tensor2_3D( hid_t group_id, const char * data_name,
        const std::array<double, 9> &value ) const;
    void write_Tensor2_3D( const HDF5_Group &group, const char * data_name,
        const std::array<double, 9> &value ) const
    { write_Tensor2_3D( group.id(), data_name, value ); }

    void write_Tensor2_3D( const char * data_name, 
        const std::array<double, 9> &value ) const;

    // --------------------------------------------------------------
    // Matrix writer
    // --------------------------------------------------------------
    void write_intMatrix( hid_t group_id,
        const char * data_name, const std::vector<int> &value,
        int row_num, int col_num ) const;
    void write_intMatrix( const HDF5_Group &group, const char * data_name,
        const std::vector<int> &value, int row_num, int col_num ) const
    { write_intMatrix( group.id(), data_name, value, row_num, col_num ); }

    void write_doubleMatrix( hid_t group_id, 
        const char * data_name, const std::vector<double> &value,
        int row_num, int col_num ) const;
    void write_doubleMatrix( const HDF5_Group &group, const char * data_name,
        const std::vector<double> &value, int row_num, int col_num ) const
    { write_doubleMatrix( group.id(), data_name, value, row_num, col_num ); }

    // --------------------------------------------------------------
    // String
    // --------------------------------------------------------------
    void write_string( const char * data_name, 
        const std::string &string_input ) const;

    void write_string( hid_t group_id, const char * data_name, 
        const std::string &string_input ) const;
    void write_string( const HDF5_Group &group, const char * data_name,
        const std::string &string_input ) const
    { write_string( group.id(), data_name, string_input ); }

  private:
    hid_t file_id;

    void check_error(herr_t status, const char * funname ) const
    {
      if(status < 0) 
      {
        std::cerr<<"Error: HDF5_Writer::"<<funname<<std::endl;
        exit( EXIT_FAILURE );
      }
    }

    void write_string_impl(hid_t location_id, const char * data_name, 
        const std::string& string_input ) const;

    void write_intScalar_impl( hid_t location_id, const char * data_name, 
        int value ) const; 

    void write_doubleScalar_impl( hid_t location_id, const char * data_name, 
        double value ) const; 

    void write_intVector_impl( hid_t location_id, const char * data_name, 
        const int * value, int length ) const;

    void write_doubleVector_impl( hid_t location_id, const char * data_name,
        const double * value, int length ) const;
};

#endif
