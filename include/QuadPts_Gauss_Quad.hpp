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

class QuadPts_Gauss_Quad final : public IQuadPts
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
   
    ~QuadPts_Gauss_Quad() override = default;

    void print_info() const override;

    int get_dim() const override {return 2;}

    int get_num_quadPts() const override {return num_pts;}

    int get_num_quadPts_x() const override {return num_pts_x;}

    int get_num_quadPts_y() const override {return num_pts_y;}

    double get_qp(const int &ii, const int &comp) const override 
    {return qp[2*ii+comp];}

    double get_qw(const int &ii) const override {return qw[ii];}

  private:
    const int num_pts, num_pts_x, num_pts_y;

    // qp : length 2 x num_pts. Stores the r-s coordinates of the
    //      quadrature points.
    // qw : length num_pts. Stores the quadrature weights.
    std::vector<double> qp {};
    std::vector<double> qw {};
};

#endif
