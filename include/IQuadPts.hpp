#ifndef IQuadPts_hpp
#define IQuadPts_hpp
// ============================================================================
// IQuadPts.hpp
// Interface class for quadrature points.
//
// Date Created: Sept. 24th 2013
// Date Modified: Jan. 17 2017
// ============================================================================
#include "Sys_Tools.hpp"

class IQuadPts
{
  public:
    IQuadPts() = default;
    
    virtual ~IQuadPts() = default;

    virtual void print_info() const
    {
      if( get_dim() == 4 )
      {
        for(int ii=0; ii<get_num_quadPts(); ++ii) SYS_T::commPrint( "%d \t qw = %e \t qp = [%e \t %e \t %e \t %e] \n", 
            ii, get_qw(ii), get_qp(ii,0), get_qp(ii,1), get_qp(ii,2), get_qp(ii,3) );
      }
      else if( get_dim() == 3 )
      {
        for(int ii=0; ii<get_num_quadPts(); ++ii) SYS_T::commPrint( "%d \t qw = %e \t qp = [%e \t %e \t %e] \n", 
            ii, get_qw(ii), get_qp(ii,0), get_qp(ii,1), get_qp(ii,2) );
      }
      else if( get_dim() == 2 )
      {
        for(int ii=0; ii<get_num_quadPts(); ++ii) SYS_T::commPrint( "%d \t qw = %e \t qp = [%e \t %e] \n", 
            ii, get_qw(ii), get_qp(ii,0), get_qp(ii,1) );
      }
      else
      {
        SYS_T::print_fatal("Error: get_dim() = %d has not been implemented.\n", get_dim() );
      }
    }
   
    // ------------------------------------------------------------------------
    // get_num_quadPts : returns the number of quadrature points 
    // ------------------------------------------------------------------------
    virtual int get_num_quadPts() const = 0;

    // ------------------------------------------------------------------------
    // get_num_quadPts_x : returns the number of quadrature points in x-direction
    //                     It is implemented for hex and quad elements.
    // ------------------------------------------------------------------------
    virtual int get_num_quadPts_x() const
    {
      SYS_T::print_fatal("Error: IQuadPts::get_num_quadPts_x is not implemented.\n");
      return -1;
    }

    // ------------------------------------------------------------------------
    // get_num_quadPts_y : returns the number of quadrature points in y-direction
    //                     It is implemented for hex and quad elements.
    // ------------------------------------------------------------------------
    virtual int get_num_quadPts_y() const
    {
      SYS_T::print_fatal("Error: IQuadPts::get_num_quadPts_y is not implemented.\n");
      return -1;
    }

    // ------------------------------------------------------------------------
    // get_num_quadPts_z : returns the number of quadrature points in z-direction
    //                     It is implemented for hex and quad elements.
    // ------------------------------------------------------------------------
    virtual int get_num_quadPts_z() const
    {
      SYS_T::print_fatal("Error: IQuadPts::get_num_quadPts_z is not implemented.\n");
      return -1;
    }

    // ------------------------------------------------------------------------
    // get_dim : returns the quadrature rule's dimension, e.g.
    //           in QuadPts, is 1,
    //           in QuadPts_Gauss_Triangle, is 3 (area coordinate),
    //           in QuadPts_Gauss_Tet, is 4 (area coordinate).
    //           qp vector will have size get_dim() x get_num_quadPts().
    // ------------------------------------------------------------------------
    virtual int get_dim() const = 0;

    // ------------------------------------------------------------------------
    // get_qp : access the ii-th quadrature point, for dim = 1 classes.    
    //          0 <= ii < get_num_quadPts()
    // ------------------------------------------------------------------------
    virtual double get_qp(const int &ii) const
    {
      SYS_T::print_fatal("Error: IQuadPts::get_qp is not implemented.\n");
      return 0.0;
    }
    
    // ------------------------------------------------------------------------
    // get_qw : access the ii-th quadrature weight
    //          0 <= ii < get_num_quadPts()
    // ------------------------------------------------------------------------
    virtual double get_qw(const int &ii) const
    {
      SYS_T::print_fatal("Error: IQuadPts::get_qw is not implemented. \n");
      return 0.0;
    }

    // ------------------------------------------------------------------------
    // get_qp : access the ii-th quadrature point's comp-th component
    //          This function is designed for quadrature points in 
    //          multidimensions.
    //          0 <= ii < get_num_quadPts()
    //          0 <= comp < get_dim()
    // ------------------------------------------------------------------------
    virtual double get_qp(const int &ii, const int &comp) const
    {
      SYS_T::print_fatal("Error: IQuadPts::get_qp is not implemented.\n");
      return 0.0;
    }

    // For User-defined QuadPts
    virtual void set_qp(const double &xi, const double &eta)
    {
      SYS_T::print_fatal("Error: IQuadPts::set_qp is not implemented.\n");
    }

    virtual void reset()
    {
      SYS_T::print_fatal("Error: IQuadPts::reset is not implemented. \n");
    }

    virtual bool check_qp_bound() const
    {
      SYS_T::print_fatal("Error: IQuadPts::check_qp_bound is not implemented. \n");
      return false;
    }
};

#endif
