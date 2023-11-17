#ifndef QUAD_PTS_BC1_HPP
#define QUAD_PTS_BC1_HPP
// ==================================================================
// QuadPts_bc1.hpp
// This is the quadrature point class that stores the sampling point
// at one end of the (0, 1) domain -- 1 point.
//
// Date: June 15 2015
// ==================================================================
#include "IQuadPts.hpp"

class QuadPts_bc1 : public IQuadPts
{
  public:
    QuadPts_bc1() = default;

    virtual ~QuadPts_bc1() = default;

    virtual void print_info() const
    {
      std::cout<<std::endl;
      std::cout<<"====== BC 1 Points ======"<<std::endl;
      std::cout<<1.0<<'\n';
      std::cout<<1.0<<'\n';
      std::cout<<"========================="<<std::endl;
    }

    virtual int get_dim() const {return 1;}

    virtual int get_num_quadPts() const {return 1;}

    virtual double get_qp(unsigned int ii) const {return 1.0;}

    virtual double get_qw(unsigned int ii) const {return 1.0;}
};

#endif
