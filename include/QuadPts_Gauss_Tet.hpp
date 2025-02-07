#ifndef QUADPTS_GAUSS_TET_HPP
#define QUADPTS_GAUSS_TET_HPP
// ==================================================================
// QuadPts_Gauss_Tet.hpp
// The Gaussian quadrature rule for a tetrahedral domain defined by 
// four vertex points:
// [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]
// 
// num_pts =  4, exact for quadratic polynomials, order of prec. 2
// num_pts =  5, exact for cubic     polynomials, order of prec. 3
// num_pts = 17, exact for quintic   polynomials, order of prec. 5
// num_pts = 29, exact for sextic    polynomials, order of prec. 6
//
// Reference: T.J.R. Hughes FEM Book p.174
//            J. Yu Symmetric Gaussian Quadrature Formulae for 
//            Tetrahedronal Regions, CMAME 43 1984:349-353
//
// Date Created: Jan. 18 2017
// ==================================================================
#include "Vec_Tools.hpp"
#include "IQuadPts.hpp"

class QuadPts_Gauss_Tet final : public IQuadPts
{
  public:
    QuadPts_Gauss_Tet( const int &in_num_pts );

    ~QuadPts_Gauss_Tet() override = default;

    void print_info() const override;

    int get_dim() const override {return 4;}

    int get_num_quadPts() const override {return num_pts;}

    double get_qp(const int &ii, const int &comp) const override
    {return qp[4*ii+comp];}

    double get_qw(const int &ii) const override
    {return qw[ii];}

  private:
    const int num_pts;

    // qp : length 4 * num_pts. Stores the r-s-t-u coordinates of the 
    //      quadrature points.
    //      u = 1 - r - s - t
    // qw : length num_pts. Stores the quadrature weights.
    std::vector<double> qp {};
    std::vector<double> qw {};
    
    // gen_permutations : generate permutations of a, b, c such that the
    //                    vector out includes the following 12 patterns.
    //                    a b c c; a c b c; a c c b;
    //                    b a c c; b c a c; b c c a;
    //                    c a b c; c b a c;
    //                    c a c b; c b c a;
    //                    c c a b; c c b a;
    std::vector<double> gen_permutations( const double &a, 
        const double &b, const double &c ) const;
};

#endif
