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

    constexpr Vector_3( const Vector_3 &source ) noexcept : vec(source.vec) {}

    constexpr Vector_3( double v0, double v1, double v2 ) 
      : vec{{ v0, v1, v2 }} {}

    constexpr Vector_3( const std::array<double,3> input ) noexcept : vec(input) {}

    ~Vector_3() = default;

    // Assignment operator
    Vector_3& operator= (const Vector_3 &source) noexcept;

    // Parenthesis operator gives access to components
    inline double& operator()(int index) noexcept {return vec[index];}

    constexpr inline const double& operator()(int index) const noexcept {return vec[index];}

    // Addition and substraction operators
    friend Vector_3 operator+( const Vector_3 &left, const Vector_3 &right ) noexcept;

    friend Vector_3 operator-( const Vector_3 &left, const Vector_3 &right ) noexcept;
    
    Vector_3& operator+=( const Vector_3 &source ) noexcept;

    Vector_3& operator-=( const Vector_3 &source ) noexcept;

    Vector_3& operator*=( double val ) noexcept;

    // unary minus operator
    Vector_3 operator-() const noexcept;

    std::vector<double> to_std_vector() const 
    {return std::vector<double>(std::begin(vec), std::end(vec));}

    constexpr std::array<double,3> to_std_array() const noexcept {return vec;}

    constexpr inline const double& x() const noexcept {return vec[0];}
    inline double& x() noexcept {return vec[0];}

    constexpr inline const double& y() const noexcept {return vec[1];}
    inline double& y() noexcept {return vec[1];}

    constexpr inline const double& z() const noexcept {return vec[2];}
    inline double& z() noexcept {return vec[2];}

    void print(std::ostream& os = std::cout, const std::string& delimiter = "\t") const;

    void gen_rand(double left =-1.0, double right = 1.0);

    constexpr inline double sum() const noexcept {return vec[0]+vec[1]+vec[2];}

    inline double norm2() const noexcept {return std::sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);}
    
    // rescale vec to be norm one and return its length
    double normalize() noexcept;

    // calculate the dot-product with a given vector
    double dot_product( const Vector_3 &source ) const noexcept;

    // Return the location of the component with the largest absolute
    // value, in [0,1,2]. Return 0 if all three components are equal.
    int get_dominant_comp() const noexcept;

  private:
    std::array<double,3> vec;
};

// calculate a scalar product of a input vector
Vector_3 operator*( double val, const Vector_3 &source ) noexcept;

namespace Vec3
{
  constexpr inline Vector_3 gen_e1() noexcept { return Vector_3(1.0, 0.0, 0.0); }

  constexpr inline Vector_3 gen_e2() noexcept { return Vector_3(0.0, 1.0, 0.0); }
    
  constexpr inline Vector_3 gen_e3() noexcept { return Vector_3(0.0, 0.0, 1.0); }
  
  constexpr inline Vector_3 gen_zero() noexcept { return Vector_3(0.0, 0.0, 0.0); }

  constexpr inline Vector_3 gen_val(double val) noexcept { return Vector_3(val, val, val); }

  inline Vector_3 gen_rand( double left =-1.0, double right = 1.0 )
  {
    std::random_device rd;
    std::mt19937_64 gen( rd() );
    std::uniform_real_distribution<double> dis(left, right);
    return Vector_3( dis(gen), dis(gen), dis(gen) );
  }

  // calculate the distance between two vector by L2 norm
  double dist( const Vector_3 &a, const Vector_3 &b ) noexcept;

  // calculate the dot product of two vectors
  double dot_product( const Vector_3 &a, const Vector_3 &b ) noexcept;

  // calculate the cross product of two vectors
  Vector_3 cross_product( const Vector_3 &a, const Vector_3 &b ) noexcept;

  // return the normalized input vector, and normalization is made by its L2 norm
  Vector_3 normalize( const Vector_3 &val ) noexcept;
}

#endif
