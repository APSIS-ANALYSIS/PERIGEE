#ifndef EXPLICIT_RK_HEUN_HPP
#define EXPLICIT_RK_HEUN_HPP

#include "ITimeMethod_RungeKutta.hpp"

class ExplicitRK_HeunRK22 : public ITimeMethod_RungeKutta 
{
  public:
    ExplicitRK_HeunRK22() : ITimeMethod_RungeKutta(2, 2, false) 
    {
      setCoefficients();
    }

  protected:
    void setCoefficients() override
    {
      // c_i values
      cc[0] = 0.0;
      cc[1] = 1.0; 

      // b_i values
      bb[0] = 0.5; 
      bb[1] = 0.5;

      // a_ij values (Butcher tableau)
      aa[0 * ss + 0] = 0.0;  aa[0 * ss + 1] = 0.0;
      aa[1 * ss + 0] = 1.0;  aa[1 * ss + 1] = 0.0;
    }
};

#endif // EXPLICIT_RK_HEUN_HPP
