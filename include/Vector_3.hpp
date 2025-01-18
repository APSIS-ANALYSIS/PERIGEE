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
#include <iomanip>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random>

class Vector_3
{
  public:
    // Default constructor generates a zero vector
    constexpr Vector_3() : vec{{ 0.0, 0.0, 0.0}} {}

    Vector_3( const Vector_3 &source ) : vec(source.vec) {}

    constexpr Vector_3( const double &v0, const double &v1, const double &v2 ) 
      : vec{{ v0, v1, v2 }} {}

    constexpr Vector_3( const std::array<double,3> input ) : vec(input) {}

    ~Vector_3() = default;

    // Assignment operator
    Vector_3& operator= (const Vector_3 &source);

    // Parenthesis operator gives access to components
    inline double& operator()(const int &index) {return vec[index];}

    inline const double& operator()(const int &index) const {return vec[index];}

    // Addition and substraction operators
    friend Vector_3 operator+( const Vector_3 &left, const Vector_3 &right );

    friend Vector_3 operator-( const Vector_3 &left, const Vector_3 &right );
    
    Vector_3& operator+=( const Vector_3 &source );

    Vector_3& operator-=( const Vector_3 &source );

    Vector_3& operator*=( const double &val );

    // unary minus operator
    Vector_3 operator-() const;

    std::vector<double> to_std_vector() const 
    {return std::vector<double>(std::begin(vec), std::end(vec));}

    std::array<double,3> to_std_array() const {return vec;}

    inline const double& x() const {return vec[0];}
    inline double& x() {return vec[0];}

    inline const double& y() const {return vec[1];}
    inline double& y() {return vec[1];}

    inline const double& z() const {return vec[2];}
    inline double& z() {return vec[2];}

    void print(std::ostream& os = std::cout, const std::string& delimiter = "\t") const;

    void gen_rand(const double &left =-1.0, const double &right = 1.0);

    inline double sum() const {return vec[0]+vec[1]+vec[2];}

    inline double norm2() const {return std::sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);}
    
    // rescale vec to be norm one and return its length
    double normalize();

    // calculate the dot-product with a given vector
    double dot_product( const Vector_3 &source ) const;

    // Return the location of the component with the largest absolute
    // value, in [0,1,2]. Return 0 if all three components are equal.
    int get_dominant_comp() const;

  private:
    std::array<double,3> vec;
};

// calculate a scalar product of a input vector
Vector_3 operator*( const double &val, const Vector_3 &source );

namespace Vec3
{
  inline Vector_3 gen_e1() { return Vector_3(1.0, 0.0, 0.0); }

  inline Vector_3 gen_e2() { return Vector_3(0.0, 1.0, 0.0); }
    
  inline Vector_3 gen_e3() { return Vector_3(0.0, 0.0, 1.0); }
  
  inline Vector_3 gen_zero() { return Vector_3(0.0, 0.0, 0.0); }

  inline Vector_3 gen_val(const double &val) { return Vector_3(val, val, val); }

  inline Vector_3 gen_rand( const double &left =-1.0, const double &right = 1.0 )
  {
    std::random_device rd;
    std::mt19937_64 gen( rd() );
    std::uniform_real_distribution<double> dis(left, right);
    return Vector_3( dis(gen), dis(gen), dis(gen) );
  }

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
