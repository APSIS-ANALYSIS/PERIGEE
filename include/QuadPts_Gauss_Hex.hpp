#ifndef QUADPTS_GAUSS_HEX_HPP
#define QUADPTS_GAUSS_HEX_HPP
// ==================================================================
// QuadPts_Gauss_Hex.hpp
// The Gaussian quadrature rule for a Hexagon domain defined by 
// eight vertex points:
// [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
// [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
//
// Date Created: Sep 7 2023
// ==================================================================
#include "Vec_Tools.hpp"
#include "IQuadPts.hpp"

class QuadPts_Gauss_Hex : public IQuadPts
{
  public:
    QuadPts_Gauss_Hex( const int &in_num_pts );

    virtual ~QuadPts_Gauss_Hex();

    virtual void print_info() const;

    virtual int get_dim() const {return 3;}

    virtual int get_num_quadPts() const {return num_pts;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const
    {return qp[3*ii+comp];}

    virtual double get_qw(unsigned int ii) const
    {return qw[ii];}

  private:
    const int num_pts;

    // qp : length 3 * num_pts. Stores the r-s-t coordinates of the 
    //      quadrature points.
    // qw : length num_pts. Stores the quadrature weights.
    std::vector<double> qp, qw;
    
    // gen_permutations : generate permutations of a, b, c such that the
    //                    vector out includes the following 12 patterns.
    //                    a b c c; a c b c; a c c b;
    //                    b a c c; b c a c; b c c a;
    //                    c a b c; c b a c;
    //                    c a c b; c b c a;
    //                    c c a b; c c b a;
    void gen_permutations( const double &a, const double &b, 
        const double &c, std::vector<double> &out ) const;
};

#endif
