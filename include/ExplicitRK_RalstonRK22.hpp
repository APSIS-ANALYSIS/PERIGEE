#ifndef EXPLICIT_RK_RALSTON_HPP
#define EXPLICIT_RK_RALSTON_HPP

#include "ITimeMethod_RungeKutta.hpp"

class ExplicitRK_RalstonRK22 : public ITimeMethod_RungeKutta 
{
  public:
    ExplicitRK_RalstonRK22() : ITimeMethod_RungeKutta(2, 2, false) 
    {
      setCoefficients();
    }

  protected:
    void setCoefficients() override
    {
      // c_i values
      cc[0] = 0.0;
      cc[1] = 2.0 / 3.0;

      // b_i values
      bb[0] = 1.0 / 4.0;
      bb[1] = 3.0 / 4.0; 

      // a_ij values (Butcher tableau)
      aa[0 * ss + 0] = 0.0;     aa[0 * ss + 1] = 0.0;
      aa[1 * ss + 0] = 2.0 / 3.0;  aa[1 * ss + 1] = 0.0;
    }
};

#endif // EXPLICIT_RK_RALSTON_HPP
