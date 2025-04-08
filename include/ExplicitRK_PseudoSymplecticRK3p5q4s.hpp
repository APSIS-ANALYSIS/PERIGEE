#ifndef EXPLICITRK_PSEUDOSYMPLECTICRK3P5Q4S_HPP
#define EXPLICITRK_PSEUDOSYMPLECTICRK3P5Q4S_HPP
// ==================================================================
// Ref: F. Capuano, G. Coppola, L. Rández, L. de Luca. Explicit Runge–Kutta schemes for 
// incompressible flow with improved energy-conservation properties. JCP 2017
// ==================================================================
#include "ITimeMethod_RungeKutta.hpp"

class ExplicitRK_PseudoSymplecticRK3p5q4s : public ITimeMethod_RungeKutta 
{
  public:
  ExplicitRK_PseudoSymplecticRK3p5q4s() : ITimeMethod_RungeKutta(4, 3, false) 
    {
      setCoefficients();
    }

  protected:
    void setCoefficients() override
    {
      const double c3 = 1.0/4.0;
   
      // c_i values 
      cc[0] = 0.0;
      cc[1] = (c3 - 1.0) / (4.0 * c3 - 3.0);
      cc[2] = c3;
      cc[3] = 1.0;

      // b_i values
      bb[0] = 1.0 / (12.0 * (c3 - 1.0));
      bb[1] = ((4.0 * c3 - 3.0) * (4.0 * c3 - 3.0)) / (12.0 * (c3 - 1.0) * (2.0 * c3 - 1.0));
      bb[2] = -1.0 / (12.0 * (c3 - 1.0) * (2.0 * c3 - 1.0));
      bb[3] = (4.0 * c3 - 3.0) / (12.0 * (c3 - 1.0));

      // a_ij values (Butcher tableau)
      aa[0 * ss + 0] = 0.0;
      aa[1 * ss + 0] = (c3 - 1.0) / (4.0 * c3 - 3.0);
      aa[2 * ss + 0] = c3 - ((2.0 * c3 - 1.0) * (4.0 * c3 - 3.0)) / (2.0 * (c3 - 1.0)); 
      aa[2 * ss + 1] = ((2.0 * c3 - 1.0) * (4.0 * c3 - 3.0)) / (2.0 * (c3 - 1.0));
      aa[3 * ss + 0] = - ((2.0 * c3 - 1.0) * (2.0 * c3 - 1.0)) / (2.0 * (c3 - 1.0)*(4.0 * c3 - 3.0)); 
      aa[3 * ss + 1] = (6.0 * c3 * c3 - 8.0 * c3 + 3.0) / (2.0 * (c3 - 1.0) * (2.0 * c3 - 1.0));
      aa[3 * ss + 2] = (c3 - 1.0) / ((2.0 * c3 - 1.0) * (4.0 * c3 - 3.0));
    }
};

#endif
