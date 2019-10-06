#ifndef HDF5_READER_HPP
#define HDF5_READER_HPP
// ==================================================================
// HDF5_Reader.hpp
// ------------------------------------------------------------------
// This is a set of tools that read in files in .h5 format. This will 
// ease the input of data from the disk.
//
// Date: July 2 2015
// ==================================================================
#include "Vec_Tools.hpp"

class HDF5_Reader
{
  public:
    // --------------------------------------------------------------
    // ! HDF5_Reader
    //   Constructor. Assuming the fild_id has been created. The 
    //   constructor will pass the given fild_id in to make it its 
    //   own variable.
    // --------------------------------------------------------------
    HDF5_Reader( const hid_t &in_file_id );
    
    
    // --------------------------------------------------------------
    // ! ~HDF5_Reader : Destructor.
    // --------------------------------------------------------------
    virtual ~HDF5_Reader();

    // --------------------------------------------------------------
    // !check_data: return a bool value that determines if the data
    //              with the specified name exists in the file.
    // --------------------------------------------------------------
    bool check_data( const char * const &name ) const;
    
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
    void read_intVector( const char * const &group_name, 
        const char * const &data_name,
        std::vector<int> &out ) const;

    
    // --------------------------------------------------------------
    // ! read_doubleVector: output the 1D integer array data into 
    //                      vector<double>.
    // --------------------------------------------------------------
    void read_doubleVector( const char * const &group_name, 
        const char * const &data_name,
        std::vector<double> &out ) const;


    // --------------------------------------------------------------
    // ! read_intMatrix : output a 2D integer Matrix into vector<int>,
    //                    with size num_row x num_col. The matrix is
    //                    stroed by rows.
    // --------------------------------------------------------------
    void read_intMatrix( const char * const &group_name,
        const char * const &data_name,
        std::vector<int> &out, int &num_row, int &num_col ) const;


    // --------------------------------------------------------------
    // ! read_doubleMatrix : output a 2D double Matrix into 
    //                       vector<double>, with size num_row x num_col. 
    //                       The matrix is stroed by rows.
    // --------------------------------------------------------------
    void read_doubleMatrix( const char * const &group_name,
        const char * const &data_name,
        std::vector<double> &out, int &num_row, int &num_col ) const;


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
    //   \para data_name : the name of the dataset in the group
    //   \para string_out : the string data
    // --------------------------------------------------------------
    void read_string( const char * const &group_name,
        const char * const &data_name,
        std::string &string_out ) const;


  private:
    const hid_t file_id;

    void check_error(const herr_t &status, const char * const &funname ) const
    {
      if( status < 0 )
      {
        std::ostringstream ss;
        ss<<"Error: HDF5_Reader::"<<funname
          <<" : status ="<<status<<" !"<<std::endl;
        SYS_T::print_exit( ss.str().c_str() );
      }
    }
};

#endif
