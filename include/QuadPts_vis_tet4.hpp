#ifndef QUADPTS_VIS_TET4_HPP
#define QUADPTS_VIS_TET4_HPP
// ==================================================================
// QuadPts_vis_tet4.hpp
//
// This is a class that stores the visualization sampling points in
// a reference tetrahedron.
//
// We use four points at [0,0,0], [1,0,0], [0,1,0], [0,0,1].
// They are the vertex points of the tetrahedron.
// 
// Note: We store them in area-coordinates like what we did in
//       QuadPts_Gauss_Tet class, so the dim = 4.
//
// Date Created: Jan. 25 2017
// ==================================================================
#include "IQuadPts.hpp"

class QuadPts_vis_tet4 : public IQuadPts
{
  public:
    QuadPts_vis_tet4();

    virtual ~QuadPts_vis_tet4() = default;

    virtual void print_info() const;

    virtual int get_dim() const {return 4;}

    virtual int get_num_quadPts() const {return 4;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const 
    {return qp[4*ii+comp];}

    virtual double get_qw(unsigned int ii) const {return qw[ii];}

  private:
    double qp[16];
    double qw[4];
};


#endif
