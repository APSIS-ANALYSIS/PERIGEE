#ifndef IQuadPts_hpp
#define IQuadPts_hpp
// ========================================================
// IQuadPts.hpp
// Interface class for quadrature points.
//
// Date Created: Sept. 24th 2013
// Date Modified: Jan. 17 2017
// ========================================================
#include "Sys_Tools.hpp"

class IQuadPts
{
  public:
    IQuadPts() = default;
    
    virtual ~IQuadPts() = default;

    virtual void print_info() const = 0;
   
    // get_num_quadPts : returns the number of quadrature points 
    virtual int get_num_quadPts() const = 0;

    // get_dim : returns the quadrature rule's dimension, e.g.
    //           in QuadPts, is 1,
    //           in QuadPts_Gauss_Triangle, is 3 (area coordinate),
    //           in QuadPts_Gauss_Tet, is 4 (area coordinate).
    //           qp vector will have size get_dim() x get_num_quadPts().
    virtual int get_dim() const = 0;

    // get_qp : access the ii-th quadrature point, for dim = 1 classes.    
    //          0 <= ii < get_num_quadPts()
    virtual double get_qp(unsigned int ii) const
    {
      SYS_T::print_fatal("Error: IQuadPts::get_qp is not implemented.\n");
      return 0.0;
    }
    
    // get_qw : access the ii-th quadrature weight
    //          0 <= ii < get_num_quadPts()
    virtual double get_qw(unsigned int ii) const
    {
      SYS_T::print_fatal("Error: IQuadPts::get_qw is not implemented. \n");
      return 0.0;
    }

    // get_qp : access the ii-th quadrature point's comp-th component
    //          This function is designed for quadrature points in 
    //          multidimensions.
    //          0 <= ii < get_num_quadPts()
    //          0 <= comp < get_dim()
    virtual double get_qp(unsigned int ii, unsigned int comp) const
    {
      SYS_T::print_fatal("Error: IQuadPts::get_qp is not implemented.\n");
      return 0.0;
    }
};

#endif
