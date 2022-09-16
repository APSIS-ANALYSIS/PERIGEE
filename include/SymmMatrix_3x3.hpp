#ifndef SYMMMATRIX_3X3_HPP
#define SYMMMATRIX_3X3_HPP
// ============================================================================
// SymmMatrix_3x3.hpp
// This is a 3 x 3 Symmetric Matrix class. The components are stored in a 1-D
// array following the Voigt notation. The matrix is thus
// 
//                    mat[0], mat[5], mat[4]
//                    mat[5], mat[1], mat[3]
//                    mat[4], mat[3], mat[2]
//
// Author: Ju Liu
// Date: Sept. 16 2022
// ============================================================================

class SymmMatrix_3x3
{
  public:
    // Constructor (default an identity 3-by-3 matrix)
    SymmMatrix_3x3();

    // Copy constructor
    SymmMatrix_3x3( const SymmMatrix_3x3 &source );

    // Constructor by six numbers in Voigt numbering
    SymmMatrix_3x3( const double &m0, const double &m1, const double &m2,
        const double &m3, const double &m4, const double &m5 );

    // Destructor
    ~SymmMatrix_3x3();

    // Parenthesis operator. It allows accessing and assigning the matrix entries.
    double& operator()(const int &index) {return mat[index];}

    const double& operator()(const int &index) const {return mat[index];}

  private:
    double mat[6];
};

#endif
