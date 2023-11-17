#ifndef QUADPTS_GAUSS_TRIANGLE_HPP
#define QUADPTS_GAUSS_TRIANGLE_HPP
// ==================================================================
// QuadPts_Gauss_Triangle.hpp
// The Gaussian quadrature rule for a triangular domain defined by vertices
// [0.0, 0.0], [1.0, 0.0], & [0.0, 1.0].
// 
// num_pts =  3, deg. of precision 2
// num_pts =  4, deg. of precision 3
// num_pts =  6, deg. of precision 4
// num_pts = 13, deg. of precision 7
//
// Reference: T.J.R.Hughes FEM book, pp 173-174.
// 
// num_pts = 19, deg. of precision 8
// num_pts = 37, deg. of precision 13
//
// DCUTRI: an algorithm for adaptive cubature over a collection of triangles, 
// ACM Transactions on Mathematical Software,
// Volume 18, Number 3, September 1992, pages 329-342.
//
// https://people.sc.fsu.edu/~jburkardt/datasets/
// quadrature_rules_tri/quadrature_rules_tri.html
//
// Date Created: Jan. 17 2017
// ==================================================================
#include "Vec_Tools.hpp"
#include "IQuadPts.hpp"

class QuadPts_Gauss_Triangle : public IQuadPts
{
  public:
    QuadPts_Gauss_Triangle( const int &in_num_pts );
    
    virtual ~QuadPts_Gauss_Triangle() = default;

    virtual void print_info() const;

    virtual int get_dim() const {return 3;}

    virtual int get_num_quadPts() const {return num_pts;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const
    {return qp[3*ii+comp];}

    virtual double get_qw(unsigned int ii) const
    {return qw[ii];}

  private:
    const int num_pts;

    // qp : length 3 x num_pts. Stores the r-s-t coordinates of the
    //      quadrature points.
    //      t = 1 - r - s.
    // qw : length num_pts. Stores the quadrature weights.
    std::vector<double> qp {};
    std::vector<double> qw {};
};

#endif
