#ifndef QUADPTS_GAUSS_QUAD_HPP
#define QUADPTS_GAUSS_QUAD_HPP
// ============================================================================
// QuadPts_Gauss_Quad.hpp
// The Gaussian quadrature rule for a quadrilateral domain defined by 
//                  [r_min, r_max] x [s_min, s_max]
//
// Date Created: Sep. 7 2023
// ============================================================================
#include "QuadPts_Gauss_1D.hpp"

class QuadPts_Gauss_Quad : public IQuadPts
{
  public:
    // Construct a quadrature rule with in_num_pts_x points in the r-direction 
    // and in_num_pts_y in the s-direction 
    QuadPts_Gauss_Quad( const int &in_num_pts_x, const int &in_num_pts_y, 
        const double &r_min = 0.0, const double &r_max = 1.0, 
        const double &s_min = 0.0, const double &s_max = 1.0 );
    
    // Construct a quadrature rule with given number of quadrature points in
    // both directions.
    QuadPts_Gauss_Quad( const int &in_num_pts_1d, 
        const double &r_min = 0.0, const double &r_max = 1.0, 
        const double &s_min = 0.0, const double &s_max = 1.0 );
   
    virtual ~QuadPts_Gauss_Quad() = default;

    virtual void print_info() const;

    virtual int get_dim() const {return 2;}

    virtual int get_num_quadPts() const {return num_pts;}

    virtual int get_num_quadPts_x() const {return num_pts_x;}

    virtual int get_num_quadPts_y() const {return num_pts_y;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const
    {return qp[2*ii+comp];}

    virtual double get_qw(unsigned int ii) const {return qw[ii];}

  private:
    const int num_pts, num_pts_x, num_pts_y;

    // qp : length 2 x num_pts. Stores the r-s coordinates of the
    //      quadrature points.
    // qw : length num_pts. Stores the quadrature weights.
    std::vector<double> qp {};
    std::vector<double> qw {};
};

#endif
