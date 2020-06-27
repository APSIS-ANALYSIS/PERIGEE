#ifndef VECTOR_3_HPP
#define VECTOR_3_HPP
// ==================================================================
// Vector_3.hpp
//
// This is a 3-component vector clas. The components are stored in
// array: vec[3]:
//                vec[0]
//                vec[1]
//                vec[2]
//
// Author: Ju Liu
// Date: June 24 2020
// ==================================================================
#include <ctime>
#include "Math_Tools.hpp"

class Vector_3
{
  public:
    // Default constructor generates a zero vector
    Vector_3();

    Vector_3( const Vector_3 &source );

    Vector_3( const double &v0, const double &v1, const double &v2 );

    ~Vector_3();

    // Copy
    void copy( const Vector_3 &source );

    void copy( double source[3] );

    // Assignment operator
    Vector_3& operator= (const Vector_3 &source);

    // Parenthesis operator gives access to components
    double& operator()(const int &index) {return vec[index];}

    const double& operator()(const int &index) const {return vec[index];}

    // Addition and substraction operators
    friend Vector_3 operator+( const Vector_3 &left, const Vector_3 &right);

    friend Vector_3 operator-( const Vector_3 &left, const Vector_3 &right);

    Vector_3& operator+=( const Vector_3 &source );

    Vector_3& operator-=( const Vector_3 &source );

    Vector_3& operator*=( const double &val );

    void print() const;

    void gen_zero();

    void gen_val(const double &val);

    void gen_rand();

    void scale( const double &val );

    double sum() const {return vec[0]+vec[1]+vec[2];}

    double norm2() const {return MATH_T::norm2(vec[0],vec[1],vec[2]);}

    void normalize(); // rescale vec to be norm one

    double dot_product( const Vector_3 &source ) const;

    friend double dot_product( const Vector_3 &a, const Vector_3 &b );

    friend Vector_3 cross_product( const Vector_3 &a, const Vector_3 &b );

  private:
    double vec[3];
};

#endif
