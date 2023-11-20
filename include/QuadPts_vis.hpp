#ifndef QUAD_PTS_VIS_HPP
#define QUAD_PTS_VIS_HPP
// ==================================================================
// QuadPts_vis.hpp
// This is the class that stores the visualization sampling points,
// which evenly distribute in (0, 1), and we require the points has
// to be more than 2, i.e., sample this element at 0 and 1.
//
// Date: Dec 18 2013
// ==================================================================
#include "IQuadPts.hpp"
#include "Vec_Tools.hpp"

class QuadPts_vis : public IQuadPts
{
  public:
    QuadPts_vis( const int &input_num_pts );
    virtual ~QuadPts_vis() = default;

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
