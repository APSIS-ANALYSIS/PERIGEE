#ifndef QUADPTS_GAUSS_HEX_HPP
#define QUADPTS_GAUSS_HEX_HPP
// ==================================================================
// QuadPts_Gauss_Hex.hpp
// The Gaussian quadrature rule for a Hexagon domain defined by 
//         [r_min, r_max] x [s_min, s_max] x [t_min, t_max]
//
// Date Created: Sep 7 2023
// ==================================================================
#include "QuadPts_Gauss_1D.hpp"

class QuadPts_Gauss_Hex final : public IQuadPts
{
  public:
    // Construct a quadrature rule with in_num_pts_x points in the r-direction, 
    // in_num_pts_y in the s-direction, and in_num_pts_z in the t-direction
    QuadPts_Gauss_Hex( const int &in_num_pts_x, 
        const int &in_num_pts_y, const int &in_num_pts_z,
        const double &r_min = 0.0, const double &r_max = 1.0,
        const double &s_min = 0.0, const double &s_max = 1.0,
        const double &t_min = 0.0, const double &t_max = 1.0 );
    
    // Construct a quadrature rule with given number of quadrature points in
    // all directions.
    QuadPts_Gauss_Hex( const int &in_num_pts_1d, 
        const double &r_min = 0.0, const double &r_max = 1.0, 
        const double &s_min = 0.0, const double &s_max = 1.0,
        const double &t_min = 0.0, const double &t_max = 1.0 );

    ~QuadPts_Gauss_Hex() override = default;

    void print_info() const override;

    int get_dim() const override {return 3;}

    int get_num_quadPts() const override {return num_pts;}

    int get_num_quadPts_x() const override {return num_pts_x;}

    int get_num_quadPts_y() const override {return num_pts_y;}

    int get_num_quadPts_z() const override {return num_pts_z;}

    double get_qp(const int &ii, const int &comp) const override
    {return qp[3*ii+comp];}

    double get_qw(const int &ii) const override
    {return qw[ii];}

  private:
    const int num_pts, num_pts_x, num_pts_y, num_pts_z;

    // qp : length 3 * num_pts. Stores the r-s-t coordinates of the 
    //      quadrature points.
    // qw : length num_pts. Stores the quadrature weights.
    std::vector<double> qp {};
    std::vector<double> qw {};
};

#endif
