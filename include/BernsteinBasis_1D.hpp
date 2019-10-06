#ifndef BERNSTEINBASIS_1D_HPP
#define BERNSTEINBASIS_1D_HPP
// ============================================================================
// BernsteinBasis_1D.hpp
//
// Bernstein Basis function defined on 1D domain [0, 1].
//  
// b_i,p = C^p_i x^i (1-x)^(p-i), i = 0, ..., p.
// C^p_i is the binomial coefficient, i.e., p! / (i!(p-i)!).
//
// We calculate the p+1 polynomials at nqua points with derivatives up to the
// second order. The container is a 1-D array with length (p+1) * nqua * 3.
// der0 : 0 - (p+1)*nqua -1
// der1 : (p+1)*nqua - 2*(p+1)*nqua - 1
// der2 : 2*(p+1)*nqua - 3*(p+1)*nqua - 1
//
// Date: Aug 25 2015
// ============================================================================
#include "IBernsteinBasis.hpp"
#include "Vec_Tools.hpp"
#include "IQuadPts.hpp"

class BernsteinBasis_1D : public IBernsteinBasis
{
  public:
    BernsteinBasis_1D( const int in_deg, const IQuadPts * const &in_quaPt );

    virtual ~BernsteinBasis_1D();

    virtual int get_deg_s() const {return deg;}

    virtual int get_nQuaPts_s() const {return nqua;}

    virtual double get_B(const int &ii, const int &qua) const
    {return val[ii + qua * deg_p1];}

    virtual double get_dB_ds(const int &ii, const int &qua) const
    {return val[ii + qua * deg_p1 + vseg];}

    virtual double get_d2B_dss(const int &ii, const int &qua) const
    {return val[ii + qua * deg_p1 + 2*vseg];}

    virtual void get_B(const int &qua, std::vector<double> &vec) const;

    virtual void get_dB_ds(const int &qua, std::vector<double> &vec) const;

    virtual void get_d2B_dss(const int &qua, std::vector<double> &vec) const;
  
  private:
    const int deg;
    const int deg_p1; // deg + 1
    const int nqua;
    const int vseg; // (deg+1)*nqua

    double * val;
    
    void eval_val( const int &quaindex, const double &x );

    BernsteinBasis_1D();
};

#endif
