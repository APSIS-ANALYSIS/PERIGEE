#ifndef MATRIX_DOUBLE_HPP
#define MATRIX_DOUBLE_HPP
// ==================================================================
// Matrix_double.hpp
// Description:
// This is a general template for double type, m by n, matrix.
//
// Date: 
// Oct. 24 2013
// ==================================================================

#include <iostream>
#include <vector>
#include <cstdlib>

#include "Vec_Tools.hpp"

using namespace std;

class Matrix_double
{
  public:
    // input: number of rows and number of collumns. generate a zero 
    //        row by col matrix
    Matrix_double( const int &row_num, const int &col_num );
    
    // input: row number and collumn number, the vector stores the rows of the
    //        matrix, with size row_num * col_num
    Matrix_double( const int &row_num, const int &col_num, 
        const vector<double> &input );

    virtual ~Matrix_double();
    
    // print matrix on screen
    virtual void print() const;

    // extract row / col vector
    virtual void get_row( const int ii, vector<double> &row_vec ) const;
    virtual void get_col( const int ii, vector<double> &col_vec ) const;

    // matrix-vector multplication
    virtual void get_MbyV( const vector<double> &in_vec, 
        vector<double> &out_vec ) const; 

    virtual int get_row_num() const
    {return row_num;}

    virtual int get_col_num() const
    {return col_num;}

    virtual double get_entry(const int &ii, const int &jj) const
    {return mat[ii*col_num + jj];}

    virtual void set_entry(const int &ii, const int &jj, double val)
    { mat[ii*col_num+jj] = val; }

    virtual void Transpose(){ cout<<"Not implemented.\n";}
    
    virtual double Det() const;
    
    virtual void Invert(){ cout<<"Not Implemented. \n"; }

  protected:
    vector<double> mat;
    int row_num;
    int col_num;
};
#endif
