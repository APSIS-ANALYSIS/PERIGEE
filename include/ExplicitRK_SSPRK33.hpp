#ifndef EXPLICIT_RK_SSPRK33_HPP
#define EXPLICIT_RK_SSPRK33_HPP

#include "ITimeMethod_RungeKutta.hpp"

class ExplicitRK_SSPRK33 : public ITimeMethod_RungeKutta 
{
  public:
    ExplicitRK_SSPRK33() : ITimeMethod_RungeKutta(3, 3, false) 
    {
        setCoefficients();
    }

  protected:
    void setCoefficients() override
    {
        // c_i values
        cc[0] = 0.0;
        cc[1] = 1.0;
        cc[2] = 0.5;

        // b_i values (final combination of k_i)
        bb[0] = 1.0 / 6.0;
        bb[1] = 1.0 / 6.0;
        bb[2] = 2.0 / 3.0;

        // a_ij values (Butcher tableau)
        aa[0 * ss + 0] = 0.0;  aa[0 * ss + 1] = 0.0;  aa[0 * ss + 2] = 0.0;
        aa[1 * ss + 0] = 1.0;  aa[1 * ss + 1] = 0.0;  aa[1 * ss + 2] = 0.0;
        aa[2 * ss + 0] = 0.25; aa[2 * ss + 1] = 0.25; aa[2 * ss + 2] = 0.0;
    }
};

#endif // EXPLICIT_RK_SSPRK33_HPP
