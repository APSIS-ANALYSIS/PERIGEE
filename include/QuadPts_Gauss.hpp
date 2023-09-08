#ifndef QUADPTS_GAUSS_HPP
#define QUADPTS_GAUSS_HPP
// ============================================================================
// QuadPts_Gauss.hpp
// This is the class that gives Gauss quadrature rule for the integration on a 
// domain [left, right] with any number of points.
//
// Date: Sept. 24th 2013
// ============================================================================
#include <limits>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "IQuadPts.hpp"

class QuadPts_Gauss : public IQuadPts
{
  public:
    QuadPts_Gauss( const int &in_num_pts );
    
    virtual ~QuadPts_Gauss();

    virtual void print_info() const;

    virtual int get_dim() const {return 1;}
    
    virtual int get_num_quadPts() const {return num_pts;}
    
    virtual double get_qp(unsigned int ii) const {return qp[ii];}
   
    virtual double get_qw(unsigned int ii) const {return qw[ii];}

  private:
    // number of quadrature points
    const int num_pts;
    
    // quadrature points and weights
    std::vector<double> qp {};
    std::vector<double> qw {};
    
    // use Newton-Raphson iteration to find the Gauss quadrature
    // points-weights. This algorithm is obtained from the dealii
    // code, quadrature_lib.cc file.
    virtual void compute_npts();

};

#endif
