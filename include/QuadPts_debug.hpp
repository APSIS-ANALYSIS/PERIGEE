#ifndef QUAD_PTS_DEBUG_HPP
#define QUAD_PTS_DEBUG_HPP
// ============================================================================
// QuadPts_debug.hpp
//
// This is the quadrature point class that we set the point location by
// constructor. This class is designed to ease the debug/test code.
//
// Date: Aug. 27 2015
// ============================================================================
#include "IQuadPts.hpp"

class QuadPts_debug final : public IQuadPts
{
  public:
    QuadPts_debug( const std::vector<double> &in_qp,
        const std::vector<double> &in_qw, const int &in_dim = 1 );

    ~QuadPts_debug() override = default;

    void print_info() const override;

    int get_dim() const override {return dim;}

    int get_num_quadPts() const override {return num_pts;}

    double get_qp(const int &ii) const override {return qp[ii];}

    double get_qp(const int &ii, const int &comp) const override
    {return qp[dim*ii + comp];}

    double get_qw(const int &ii) const override {return qw[ii];}

  private:
    const int num_pts, dim;

    const std::vector<double> qp, qw;
};

#endif
