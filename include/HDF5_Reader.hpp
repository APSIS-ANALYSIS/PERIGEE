#ifndef HDF5_READER_HPP
#define HDF5_READER_HPP
// ==================================================================
// HDF5_Reader.hpp
// ------------------------------------------------------------------
// This is a set of tools that read in files in .h5 format. This will 
// ease the input of data from the disk.
//
// To construct this class, one should first obtain a hid_t data by
// calling H5Fopen. After the data are all read, one should call
// H5Fclose.
//
// A typical way of using this class is:
//
// hid_t file_id = H5Fopen(name_of_h5_file, mode, H5P_DEFAULT)
//
// HDF5_Reader * h5r = new HDF5_Reader( file_id );
//
// call read functions
//
// delete h5r;
// H5Fclose(file_id);
//
// There are two modes for H5Fopen
// H5F_ACC_RDONLY: the application will read only the file
// H5F_ACC_RDWR: the application will read and write the file.
//
// Basically, we only use H5F_ACC_RDONLY.
// For more details, see 
//      https://support.hdfgroup.org/HDF5/Tutor/crtfile.html
// 
// Author: Ju Liu
// Date: July 2 2015
// ==================================================================
#include <array>
#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "hdf5.h"

class HDF5_Reader
{
  public:
    // --------------------------------------------------------------
    // ! HDF5_Reader
    //   Constructor. Assuming the fild_id has been created. The 
    //   constructor will pass the given fild_id in to make it its 
    //   own variable.
    // --------------------------------------------------------------
    HDF5_Reader( const hid_t &in_file_id ) : file_id(in_file_id) {}
    
    // --------------------------------------------------------------
    // ! ~HDF5_Reader : Destructor.
    // --------------------------------------------------------------
    virtual ~HDF5_Reader() = default;

    // --------------------------------------------------------------
    // !check_data: return a bool value that determines if the data
    //              with the specified name exists in the file.
    // --------------------------------------------------------------
    bool check_data( const char * const &name ) const
    {
      return H5Lexists(file_id, name, H5P_DEFAULT);
    }
    
    // --------------------------------------------------------------
    // !read_intScalar: return the integer scalar data by specifing 
    //                  its group_name with data_name.
    // --------------------------------------------------------------
    int read_intScalar( const char * const &group_name,
        const char * const &data_name ) const;

    // --------------------------------------------------------------
    // !read_doubleScalar: return the double scalar data by specifing 
    //                  its group_name with data_name.
    // --------------------------------------------------------------
    double read_doubleScalar( const char * const &group_name,
        const char * const &data_name ) const;

    // --------------------------------------------------------------
    // ! read_intVector: output the 1D integer array data into 
    //                   vector<int>.
    // --------------------------------------------------------------
    std::vector<int> read_intVector( const char * const &group_name,
        const char * const &data_name ) const;
    
    // --------------------------------------------------------------
    // ! read_doubleVector: output the 1D integer array data into 
    //                      vector<double>.
    // --------------------------------------------------------------
    std::vector<double> read_doubleVector( const char * const &group_name,
        const char * const &data_name ) const;

    // --------------------------------------------------------------
    // ! read_Vector_3 : output the std::array<double, 3>.
    // --------------------------------------------------------------
    std::array<double, 3> read_Vector_3(const char* const& group_name,
        const char * const &data_name ) const;

    // --------------------------------------------------------------
    // ! read_Tensor2_3D : output the std::array<double, 9>.
    // --------------------------------------------------------------
    std::array<double, 9> read_Tensor2_3D(const char* const& group_name,
        const char * const &data_name ) const;

    // --------------------------------------------------------------
    // ! read_intMatrix : output a 2D integer Matrix into vector<int>,
    //                    with size num_row x num_col. The matrix is
    //                    stroed by rows.
    // --------------------------------------------------------------
    std::vector<int> read_intMatrix( const char * const &group_name,
        const char * const &data_name, int &num_row, int &num_col ) const;

    // --------------------------------------------------------------
    // ! read_doubleMatrix : output a 2D double Matrix into 
    //                       vector<double>, with size num_row x num_col. 
    //                       The matrix is stroed by rows.
    // --------------------------------------------------------------
    std::vector<double> read_doubleMatrix( const char * const &group_name,
        const char * const &data_name, int &num_row, int &num_col ) const;

    // --------------------------------------------------------------
    // ! read_intArray
    //   This is the read function that will read from the file 
    //   specified by file_id.
    //   \para group_name : the name of the group in the file of file_id
    //   \para data_name : the name of the dateset in the group of 
    //                     group_name
    //   \para data_rank : the dimension of the target object
    //   \para data_dims : the length in each dimesion. e.g. 
    //                     for a 3 x 5 matrix,
    //                     rank = 2 and data_dims [] = {3,5}. 
    //   \para data : the pointer to the array which passes the value 
    //                out. For two-dimensional arrays, the data will 
    //                be output in the format of a row array.
    //                e.g. [1,2;
    //                      3,4] 
    //                will become {1,2,3,4} in the output data.
    //    Note: The user is responsible for deleting/freeing the memory 
    //          space for data_dims and data. 
    // --------------------------------------------------------------
    void read_intArray(const char * const &group_name, 
        const char * const &data_name,
        hid_t &data_rank, hsize_t * &data_dims, int * &data  ) const;

    // --------------------------------------------------------------
    // ! read_doubleArray
    //   This is the read function that will read from the file 
    //   specified by file_id.
    //   \para group_name : the name of the group in the file of file_id
    //   \para data_name : the name of the dateset in the group of 
    //                     group_name
    //   \para data_rank : the dimension of the target object
    //   \para data_dims : the length in each dimesion. e.g. 
    //                     for a 3 x 5 matrix,
    //                     rank = 2 and data_dims [] = {3,5}. 
    //   \para data : the pointer to the array which passes the value 
    //                out. For two-dimensional arrays, the data will 
    //                be the row array.
    //                e.g. [1.0, 2.0; 
    //                      3.0, 4.0] 
    //                will become {1.0, 2.0, 3.0, 4.0} in the output data.
    //    Note: The user is responsible for deleting/freeing the memory 
    //          space for data_dims and data. 
    // --------------------------------------------------------------
    void read_doubleArray(const char * const &group_name, 
        const char * const &data_name,
        hid_t &data_rank, hsize_t * &data_dims, double * &data  ) const;

    // --------------------------------------------------------------
    // ! read_string
    //   This is the read function that will read from the file
    //   specified by file_id.
    //   \para group_name : the name of the group
    //   \para data_name  : the name of the dataset in the group
    // --------------------------------------------------------------
    std::string read_string( const char * const &group_name,
        const char * const &data_name ) const;

  private:
    const hid_t file_id;

    void check_error(const herr_t &status, const char * const &funname ) const
    {
      if( status < 0 )
      {
        std::ostringstream ss;
        ss<<"Error: HDF5_Reader::"<<funname
          <<" : status ="<<status<<" !"<<std::endl;
        SYS_T::print_fatal( ss.str().c_str() );
      }
    }
};

#endif
