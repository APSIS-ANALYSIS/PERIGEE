#ifndef BERNSTEIN_BASIS_HPP
#define BERNSTEIN_BASIS_HPP
// ==================================================================
// BernsteinBasis.hpp
// Bernstein polynomial class. It calculate the Bernstein polynomial
// function and derivative value up to 2nd order derivative.
// Date:
// Oct. 31 2013
// ==================================================================
#include <iostream>
#include <vector>
#include <cstdlib>
#include "Vec_Tools.hpp"

using namespace std;

class BernsteinBasis
{
  public:
    BernsteinBasis( const int &input_degree, const double &input_xi );
    virtual ~BernsteinBasis();
    
    virtual int get_degree() const {return degree;}
    virtual double get_xi() const {return xi;}

    virtual double get_der0(const int &fun_index) const {return der0[fun_index];}
    virtual double get_der1(const int &fun_index) const {return der1[fun_index];}
    virtual double get_der2(const int &fun_index) const {return der2[fun_index];}

    virtual void print() const;    
  private:
    // This constructor is not allowed.
    BernsteinBasis();
    
    int degree;
    double xi;

    // 0th, first, second order derivatives of all 
    // bernstein basis funcitions (p+1), at given location xi.
    vector<double> der0;
    vector<double> der1;
    vector<double> der2;

    // detailed computation of der0, der1, and der2.
    // calculate the function by switch the polynomial degree
    void funcInit_switch();
};

#endif
