#ifndef EXPLICITRK_EM_HPP
#define EXPLICITRK_EM_HPP

#include "ITimeMethod_RungeKutta.hpp"

class ExplicitRK_EMRK22 : public ITimeMethod_RungeKutta 
{
  public:
    ExplicitRK_EMRK22() : ITimeMethod_RungeKutta(2, 2, false) 
    {
      setCoefficients();
    }

  protected:
    void setCoefficients() override
    {
      // c_i values
      cc[0] = 0.0;
      cc[1] = 0.5;

      // b_i values
      bb[0] = 0.0;
      bb[1] = 1.0;

      // a_ij values (Butcher tableau)
      aa[0 * ss + 0] = 0.0;  aa[0 * ss + 1] = 0.0;
      aa[1 * ss + 0] = 0.5;  aa[1 * ss + 1] = 0.0;
    }
};

#endif
