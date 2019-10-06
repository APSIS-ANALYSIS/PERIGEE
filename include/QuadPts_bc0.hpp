#ifndef QUAD_PTS_BC0_HPP
#define QUAD_PTS_BC0_HPP
// ==================================================================
// QuadPts_bc0.hpp
// This is the quadrature point class that stores the sampling points
// at one end of the (0, 1) domain -- 0 point.
//
// Date: June 15 2015
// ==================================================================
#include "Vec_Tools.hpp"
#include "IQuadPts.hpp"

class QuadPts_bc0 : public IQuadPts
{
  public:
    QuadPts_bc0();
    virtual ~QuadPts_bc0();

    virtual void print_info() const;

    virtual int get_dim() const {return 1;}

    virtual int get_num_quadPts() const {return num_pts;}

    virtual double get_qp(unsigned int ii) const {return qp[ii];}

    virtual double get_qw(unsigned int ii) const {return qw[ii];}

  private:
    int num_pts;

    std::vector<double> qp, qw;
};

#endif
