#ifndef BERNSTEIN_BASIS_ARRAY_HPP
#define BERNSTEIN_BASIS_ARRAY_HPP
// ==================================================================
// BernsteinBasis_Array.hpp
// Bernstein polynomial class. It calculate the Bernstein polynomial
// function and derivative value up to 2nd order derivative at givne
// quadrature points.
//
// Date:
// Oct. 31 2013
// ==================================================================
#include "Math_Tools.hpp"
#include "IQuadPts.hpp"

class BernsteinBasis_Array
{
  public:
    // --------------------------------------------------------------
    // Construct a Bernstein polynomail array evaluated at a suite of
    // quadrature points given by the IQuadPts object.
    // --------------------------------------------------------------
    BernsteinBasis_Array( const int &input_degree, 
        const class IQuadPts * const &in_quaInfo );
    
    // --------------------------------------------------------------
    // Construct a Bernstein polynomial array evaluated at a given 
    // point: in_qp (the number of quadrature pt is 1).
    // --------------------------------------------------------------
    BernsteinBasis_Array( const int &input_degree, const double &in_qp );
   
    virtual ~BernsteinBasis_Array();
    
    virtual int get_degree() const {return degree;}

    virtual double get_nQuapts() const {return nQuapts;}


    virtual double get_der0(const int &fun_index, const int &quaindex) const
    {return der0[ degp1*quaindex + fun_index ];}
    

    virtual double get_der1(const int &fun_index, const int &quaindex) const
    {return der1[ degp1*quaindex + fun_index ];}
    

    virtual double get_der2(const int &fun_index, const int &quaindex) const
    {return der2[ degp1*quaindex + fun_index ];}

    virtual void print_info() const;    

  private:
    // This constructor is not allowed.
    BernsteinBasis_Array();
    
    const int degree;
    const int nQuapts;
    const int degp1;

    // 0th, first, second order derivatives of all 
    // bernstein basis funcitions (p+1), at given location xi.
    double * der0;
    double * der1;
    double * der2;

    // --------------------------------------------------------------
    // Calculation of der0, der1, and der2.
    // calculate the function by switch the polynomial degree
    // degree = 0, one constant basis
    // degree = 1, two linear basis
    // degree = 2, three quadratic basis
    // degree = 3, four cubic basis
    // Otherwise, basis undefined and return an error.
    // --------------------------------------------------------------
    void funcInit_switch(const int &quaindex, const double &xi);
};

#endif
