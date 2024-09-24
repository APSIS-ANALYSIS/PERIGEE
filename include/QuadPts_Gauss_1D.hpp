#ifndef QUADPTS_GAUSS_1D_HPP
#define QUADPTS_GAUSS_1D_HPP
// ============================================================================
// QuadPts_Gauss_1D.hpp
// This is the class that gives Gauss quadrature rule for the integration on a 
// domain [min, max] with any number of points.
//
// Date: Sept. 24th 2013
// ============================================================================
#include <limits>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "IQuadPts.hpp"

class QuadPts_Gauss_1D final : public IQuadPts
{
  public:
    // Construct in_num_pts-point rule for [min, max] domain
    QuadPts_Gauss_1D( const int &in_num_pts, const double &min = 0.0, const double &max = 1.0 );
    
    ~QuadPts_Gauss_1D() override = default;

    void print_info() const override;

    int get_dim() const override {return 1;}
    
    int get_num_quadPts() const override {return num_pts;}
    
    double get_qp(const int &ii) const override {return qp[ii];}
   
    double get_qw(const int &ii) const override {return qw[ii];}

  private:
    // number of quadrature points
    const int num_pts;
    
    // quadrature points and weights
    std::vector<double> qp {};
    std::vector<double> qw {};
    
    // use Newton-Raphson iteration to find the Gauss quadrature
    // points-weights. This algorithm is obtained from the dealii
    // code, quadrature_lib.cc file.
    void compute_npts();
};

#endif
