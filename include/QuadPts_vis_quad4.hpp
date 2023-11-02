#ifndef QUADPTS_VIS_QUAD4_HPP
#define QUADPTS_VIS_QUAD4_HPP
// ==================================================================
// QuadPts_vis_quad4.hpp
//
// This is a class of quadrature points for 4 node quad element.
//
// We use [0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0].
// These are the vertex points for the quad elements.
// ==================================================================
#include <vector>
#include "IQuadPts.hpp"

class QuadPts_vis_quad4 : public IQuadPts
{
  public:
    QuadPts_vis_quad4();

    virtual ~QuadPts_vis_quad4();

    virtual void print_info() const;

    // it stores the coordinate of the quadrature points 
    // in the sequence of r-s, so the dim is 2
    virtual int get_dim() const {return 2;}

    virtual int get_num_quadPts() const {return 2;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const
    {return qp[2*ii+comp];}

    virtual double get_qw(unsigned int ii) const
    {return qw[ii];}

  private:
    double qp[8];
    double qw[4];
};

#endif
