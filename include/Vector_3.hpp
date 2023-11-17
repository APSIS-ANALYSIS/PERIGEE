#ifndef VECTOR_3_HPP
#define VECTOR_3_HPP
// ============================================================================
// Vector_3.hpp
//
// This is a 3-component vector class. The components are stored in an array 
// vec[3]:
//             [ vec[0]; vec[1]; vec[2] ]
// 
// It is designed primarily for material constitutive routines.
//
// Author: Ju Liu
// Date: June 24 2020
// ============================================================================
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <array>
#include <random>

class Vector_3
{
  public:
    // Default constructor generates a zero vector
    Vector_3();

    Vector_3( const Vector_3 &source );

    Vector_3( const double &v0, const double &v1, const double &v2 );

    ~Vector_3() = default;

    // Copy
    void copy( const Vector_3 &source );

    void copy( const double source[3] );

    // Assignment operator
    Vector_3& operator= (const Vector_3 &source);

    // Parenthesis operator gives access to components
    double& operator()(const int &index) {return vec[index];}

    const double& operator()(const int &index) const {return vec[index];}

    // Addition and substraction operators
    friend Vector_3 operator+( const Vector_3 &left, const Vector_3 &right );

    friend Vector_3 operator-( const Vector_3 &left, const Vector_3 &right );
    
    Vector_3& operator+=( const Vector_3 &source );

    Vector_3& operator-=( const Vector_3 &source );

    Vector_3& operator*=( const double &val );

    std::vector<double> to_std_vec() const;

    std::array<double, 3> to_std_array() const;

    const double& x() const {return vec[0];}
    double& x() {return vec[0];}

    const double& y() const {return vec[1];}
    double& y() {return vec[1];}

    const double& z() const {return vec[2];}
    double& z() {return vec[2];}

    void print() const;

    void gen_zero();

    void gen_val(const double &val);

    void gen_rand(const double &left =-1.0, const double &right = 1.0);

    void gen_e1() {vec[0]=1.0; vec[1]=0.0; vec[2]=0.0;}
    
    void gen_e2() {vec[0]=0.0; vec[1]=1.0; vec[2]=0.0;}
    
    void gen_e3() {vec[0]=0.0; vec[1]=0.0; vec[2]=1.0;}

    double sum() const {return vec[0]+vec[1]+vec[2];}

    double norm2() const {return std::sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);}
    
    // rescale vec to be norm one and return its length
    double normalize();

    // calculate the dot-product with a given vector
    double dot_product( const Vector_3 &source ) const;

    // Return the location of the component with the largest absolute
    // value, in [0,1,2]. Return 0 if all three components are equal.
    int get_dominant_comp() const;

  private:
    double vec[3];
};

// calculate a scalar product of a input vector
Vector_3 operator*( const double &val, const Vector_3 &source );

namespace Vec3
{
  // calculate the distance between two vector by L2 norm
  double dist( const Vector_3 &a, const Vector_3 &b );

  // calculate the dot product of two vectors
  double dot_product( const Vector_3 &a, const Vector_3 &b );

  // calculate the cross product of two vectors
  Vector_3 cross_product( const Vector_3 &a, const Vector_3 &b );

  // return the normalized input vector, and normalization is made by its L2 norm
  Vector_3 normalize( const Vector_3 &val );
}

#endif
