#ifndef BERNSTEINBASIS_2D_HPP
#define BERNSTEINBASIS_2D_HPP
// ============================================================================
// BernsteinBasis_2D.hpp
//
// Bernstein polynomial basis functions defined for 2D domains. 
// ============================================================================
#include "BernsteinBasis_1D.hpp"

class BernsteinBasis_2D : public IBernsteinBasis
{
  public:
    BernsteinBasis_2D( const int &in_sdeg, const int &in_tdeg,
        const IQuadPts * const &in_quaPt_s,
        const IQuadPts * const &in_quaPt_t );

    virtual ~BernsteinBasis_2D();

    virtual int get_deg_s() const {return sdeg;}

    virtual int get_deg_t() const {return tdeg;}

    virtual int get_nQuaPts_s() const {return nqua_s;}
    
    virtual int get_nQuaPts_t() const {return nqua_t;}

    // Get the Bernstein elements function and derivatives up to the second
    // order. 
    // ii is the index for the basis function. 
    // 0 <= ii < (sdeg+1)(tdeg+1).
    // qua is the index for the quadrature points.
    // 0 <= qua < nqua_s * nqua_t.
    virtual double get_B(const int &ii, const int &qua) const;

    virtual double get_dB_ds(const int &ii, const int &qua) const;

    virtual double get_dB_dt(const int &ii, const int &qua) const;

    virtual double get_d2B_dss(const int &ii, const int &qua) const;

    virtual double get_d2B_dst(const int &ii, const int &qua) const;

    virtual double get_d2B_dtt(const int &ii, const int &qua) const;

    // Get the Bernstein basis as a vector from 0 to (sdeg_p1)*(tdeg_p1)
    virtual void get_B(const int &qua_s, const int &qua_t, std::vector<double> &vec) const;
    virtual void get_dB_ds(const int &qua_s, const int &qua_t, std::vector<double> &vec) const;
    virtual void get_dB_dt(const int &qua_s, const int &qua_t, std::vector<double> &vec) const;
    
    virtual void get_d2B_dss(const int &qua_s, const int &qua_t, std::vector<double> &vec) const;
    virtual void get_d2B_dst(const int &qua_s, const int &qua_t, std::vector<double> &vec) const;
    virtual void get_d2B_dtt(const int &qua_s, const int &qua_t, std::vector<double> &vec) const;

    virtual void print_info() const;

  private:
    const int sdeg;
    const int tdeg;

    const int sdeg_p1; // sdeg + 1
    
    const int tdeg_p1; // tdeg + 1

    const int sp1tp1; // (sdeg+1)(tdeg+1)

    const int nqua_s;
    const int nqua_t;

    BernsteinBasis_1D * bs;
    BernsteinBasis_1D * bt;
};

#endif
