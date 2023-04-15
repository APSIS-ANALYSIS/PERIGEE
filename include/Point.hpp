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
      assert(dim == 2);
      this -> coor[0] = x;
      this -> coor[1] = y;
    }

    Point( const T &x, const T &y, const T &z )
    {
      assert(dim == 3);
      this -> coor[0] = x;
      this -> coor[1] = y;
      this -> coor[2] = z;
    }
    
    ~Point() {};

    // read-write operator
    T& operator[](const int &ii) 
    {assert(ii>=0 && ii<dim); return coor[ii];}

    // read only access operator
    const T& operator[](const int &ii) const 
    {assert(ii>=0 && ii<dim); return coor[ii];}

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
  assert(dim > 0);
  for(int ii=0; ii<dim; ++ii) this -> coor[ii] = rhs[ii];

  return *this;
}


template<int dim, typename T> inline
Point<dim,T> &
Point<dim,T>::operator += ( const Point<dim,T> &rhs )
{
  assert(dim > 0);
  for(int ii=0; ii<dim; ++ii) this -> coor[ii] += rhs[ii];

  return *this;
}


template<int dim, typename T> inline
Point<dim,T> &
Point<dim,T>::operator -= ( const Point<dim,T> &rhs )
{
  assert(dim > 0);
  for(int ii=0; ii<dim; ++ii) this -> coor[ii] -= rhs[ii];

  return *this;
}


template<int dim, typename T> inline
Point<dim,T> &
Point<dim,T>::operator *= ( const T &val )
{
  assert(dim > 0);
  for(int ii=0; ii<dim; ++ii) this -> coor[ii] *= val;

  return *this;
}


template<int dim, typename T> inline
Point<dim,T> &
Point<dim,T>::operator /= ( const T &val )
{
  assert(dim > 0 && val != 0.0);
  for(int ii=0; ii<dim; ++ii) this -> coor[ii] /= val;

  return *this;
}

#endif
