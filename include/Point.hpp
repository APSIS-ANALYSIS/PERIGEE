#ifndef POINT_HPP
#define POINT_HPP
// ==================================================================
// Point.hpp
// 
// The Point class represent a point in n-dimensional space, and this
// class stores the point's coordinates.
//
// dim: An integer that denotes the dimension of the space in which a 
// point lies.
//
// number: The data type in which the coordinates values are to be 
// stored. This will in most cases be defaulted to double type.
//
// Date Created: Jan. 18 2017
// ==================================================================
#include <iostream>

template <int dim, typename T = double>
class Point
{
  public:
    Point() {for(int ii=0; ii<dim; ++ii) this->coor[ii] = static_cast<T>(0.0);}
    
    Point( const T &x, const T &y )
    {
      ASSERT(dim == 2, "The dimension of Point<dim,T> has to equal 2.\n" );
      this -> coor[0] = x;
      this -> coor[1] = y;
    }

    Point( const T &x, const T &y, const T &z )
    {
      ASSERT(dim == 3, "The dimension of Point<dim,T> has to equal 3.\n" );
      this -> coor[0] = x;
      this -> coor[1] = y;
      this -> coor[2] = z;
    }
    
    ~Point() {};

    // read-write operator
    T& operator[](const int &ii) 
    {ASSERT(ii>=0 && ii<dim, "The read-write access operator of the Point class error.\n" ); return coor[ii];}

    // read only access operator
    const T& operator[](const int &ii) const 
    {ASSERT(ii>=0 && ii<dim, "The read only access operator of the Point class error.\n" ); return coor[ii];}

    Point& operator = ( const Point<dim, T> &rhs );

    Point& operator += ( const Point<dim, T> &rhs );

    Point& operator -= ( const Point<dim, T> &rhs );

    Point& operator *= ( const T &val );
    
    Point& operator /= ( const T &val );

  private:
    T coor[dim];
};


template<int dim, typename T> inline
std::ostream &operator << (std::ostream &out, const Point<dim,T> &p)
{
  out<<"Point(";
  for(int ii=0; ii<dim-1; ++ii) out<<p[ii]<<' ';
  out<<p[dim-1]<<")";
  return out;
}


template<int dim, typename T> inline
Point<dim,T> &
Point<dim,T>::operator = ( const Point<dim,T> &rhs )
{
  ASSERT(dim > 0, "The operator = of the Point class error, the dimension of Point<dim,T> has to bigger than 0.\n" );
  for(int ii=0; ii<dim; ++ii) this -> coor[ii] = rhs[ii];

  return *this;
}


template<int dim, typename T> inline
Point<dim,T> &
Point<dim,T>::operator += ( const Point<dim,T> &rhs )
{
  ASSERT(dim > 0, "The operator += of the Point class error, the dimension of Point<dim,T> has to bigger than 0.\n" );
  for(int ii=0; ii<dim; ++ii) this -> coor[ii] += rhs[ii];

  return *this;
}


template<int dim, typename T> inline
Point<dim,T> &
Point<dim,T>::operator -= ( const Point<dim,T> &rhs )
{
  ASSERT(dim > 0, "The operator -= of the Point class error, the dimension of Point<dim,T> has to bigger than 0.\n" );
  for(int ii=0; ii<dim; ++ii) this -> coor[ii] -= rhs[ii];

  return *this;
}


template<int dim, typename T> inline
Point<dim,T> &
Point<dim,T>::operator *= ( const T &val )
{
  ASSERT(dim > 0, "The operator *= of the Point class error, the dimension of Point<dim,T> has to bigger than 0.\n" );
  for(int ii=0; ii<dim; ++ii) this -> coor[ii] *= val;

  return *this;
}


template<int dim, typename T> inline
Point<dim,T> &
Point<dim,T>::operator /= ( const T &val )
{
  ASSERT(dim > 0 && val != 0.0, "The operator /= of the Point class error, the dimension of Point<dim,T> has to bigger than 0 and the divisor should not be 0.0\n" );
  for(int ii=0; ii<dim; ++ii) this -> coor[ii] /= val;

  return *this;
}

#endif
