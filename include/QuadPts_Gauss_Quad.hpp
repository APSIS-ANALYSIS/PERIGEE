#ifndef QUADPTS_GAUSS_QUAD_HPP
#define QUADPTS_GAUSS_QUAD_HPP
// ==================================================================
// QuadPts_Gauss_Quad.hpp
// The Gaussian quadrature rule for a quadrelateral domain defined by 
// vertices
// [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]
//
// Date Created: Sep. 7 2023
// ==================================================================
#include "Vec_Tools.hpp"
#include "IQuadPts.hpp"

class QuadPts_Gauss_Quad : public IQuadPts
{
  public:
    QuadPts_Gauss_Quad( const int &in_num_pts );
    
    virtual ~QuadPts_Gauss_Quad();

    virtual void print_info() const;

    virtual int get_dim() const {return 2;}

    virtual int get_num_quadPts() const {return num_pts;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const
    {return qp[2*ii+comp];}

    virtual double get_qw(unsigned int ii) const {return qw[ii];}

  private:
    const int num_pts;

    // qp : length 2 x num_pts. Stores the r-s coordinates of the
    //      quadrature points.
    // qw : length num_pts. Stores the quadrature weights.
    std::vector<double> qp {};
    std::vector<double> qw {};
};

#endif
