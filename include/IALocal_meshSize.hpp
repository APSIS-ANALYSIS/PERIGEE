#ifndef IALOCAL_MESHSIZE_HPP
#define IALOCAL_MESHSIZE_HPP
// ==================================================================
// IALocal_meshSize.hpp
// Interface for analysis use local mesh size class.
//
// Date:
// Nov. 11 2013
// ==================================================================
#include <iostream>

class IALocal_meshSize
{
  public:
    IALocal_meshSize(){};
    virtual ~IALocal_meshSize(){};

    virtual void print() const {std::cout<<"print() not implemented. \n";}

    // get_meshsize
    // returns the element's size in the reference domain.
    // If the element is 1D, this function returns hx;
    // If the element is 2D, this function returns hx * hy;
    // If the element is 3D, this function returns hx * hy * hz.
    virtual double get_meshsize(const int &e) const = 0;

    virtual double get_hx(const int &e) const = 0;
    virtual double get_hy(const int &e) const {return 0.0;}
    virtual double get_hz(const int &e) const {return 0.0;}
};


#endif
