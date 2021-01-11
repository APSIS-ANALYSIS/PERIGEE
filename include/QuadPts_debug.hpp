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
#include "Vec_Tools.hpp"

class QuadPts_debug : public IQuadPts
{
  public:
    QuadPts_debug( const int &len, const std::vector<double> &in_qp,
        const std::vector<double> &in_qw );

    QuadPts_debug( const int &in_dim, const int &in_numpt, 
        const std::vector<double> &in_qp, const std::vector<double> &in_qw );

    virtual ~QuadPts_debug();

    virtual void print_info() const;

    virtual int get_dim() const {return dim;}

    virtual int get_num_quadPts() const {return num_pts;}

    virtual double get_qp(unsigned int ii) const {return qp[ii];}

    virtual double get_qw(unsigned int ii) const {return qw[ii];}

  private:
    const int num_pts, dim;

    std::vector<double> qp, qw;
};

#endif
