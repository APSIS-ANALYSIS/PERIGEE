#ifndef QUADPTS_USERDEFINED_TRIANGLE_HPP
#define QUADPTS_USERDEFINED_TRIANGLE_HPP

#include "IQuadPts.hpp"

class QuadPts_UserDefined_Triangle : public IQuadPts
{
  public:
    QuadPts_UserDefined_Triangle()
    {
      reset();
    }
    
    virtual ~QuadPts_UserDefined_Triangle() = default;

    virtual int get_dim() const {return 3;}

    virtual int get_num_quadPts() const {return 1;}

    virtual double get_qp(const int &comp) const
    {return qp[comp];}

    virtual double get_qp(const int &ii, const int &comp) const
    {return qp[3 * ii + comp];}

    virtual void set_qp(const double &xi, const double &eta)
    {
      qp[0] = xi;
      qp[1] = eta;
      qp[2] = 1.0 - qp[0] - qp[1];
    }

    virtual void reset()
    {
      qp[0] = 0.333333333333333;
      qp[1] = 0.333333333333333;
      qp[2] = 0.333333333333333;
    }

    virtual bool check_qp_bound() const
    {
      const double epsilon = 2.0e-3;

      if( qp[0]>=(0-epsilon) && qp[0]<=(1+epsilon) && qp[1]>=(0-epsilon) && qp[1]<=(1+epsilon) && qp[2]>=(0-epsilon) )
        return true;
      else
        return false;
    }

  private:
    // qp : length 3. Stores the r-s-t coordinates of the
    //      quadrature point.
    //      t = 1 - r - s.
    double qp[3];
};

#endif
