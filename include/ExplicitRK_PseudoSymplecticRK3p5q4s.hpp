#ifndef EXPLICITRK_PSEUDOSYMPLECTICRK3P5Q4S_HPP
#define EXPLICITRK_PSEUDOSYMPLECTICRK3P5Q4S_HPP
// ==================================================================
// Explicit pseudo-symplectic Runge–Kutta method: 3p5q4s
// Four-stage scheme with third-order accuracy on solution and fifth-order 
// accuracy on energy conservation
//
// Ref: F. Capuano et al. Explicit Runge-Kutta schemes for incompressible flow 
//      with improved energy-conservation properties. Journal of Computational Physics, 2017
//
// Date: Apr. 14 2025
// Author: Yujie Sun
// ==================================================================
#include "ITimeMethod_RungeKutta.hpp"

class ExplicitRK_PseudoSymplecticRK3p5q4s final : public ITimeMethod_RungeKutta 
{
  public:
    ExplicitRK_PseudoSymplecticRK3p5q4s() : ITimeMethod_RungeKutta(4) 
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

      // a_ij values
      aa[0 * ss + 0] = 0.0;
      aa[1 * ss + 0] = (c3 - 1.0) / (4.0 * c3 - 3.0);
      aa[2 * ss + 0] = c3 - ((2.0 * c3 - 1.0) * (4.0 * c3 - 3.0)) / (2.0 * (c3 - 1.0)); 
      aa[2 * ss + 1] = ((2.0 * c3 - 1.0) * (4.0 * c3 - 3.0)) / (2.0 * (c3 - 1.0));
      aa[3 * ss + 0] = - ((2.0 * c3 - 1.0) * (2.0 * c3 - 1.0)) / (2.0 * (c3 - 1.0)*(4.0 * c3 - 3.0)); 
      aa[3 * ss + 1] = (6.0 * c3 * c3 - 8.0 * c3 + 3.0) / (2.0 * (c3 - 1.0) * (2.0 * c3 - 1.0));
      aa[3 * ss + 2] = (c3 - 1.0) / ((2.0 * c3 - 1.0) * (4.0 * c3 - 3.0));
    }

    ~ExplicitRK_PseudoSymplecticRK3p5q4s() override = default;

    double get_RK_a(const int &ii, const int &jj) const override
    {
      ASSERT((ii >= 0 && ii < ss && jj >= 0 && jj < ss), "Error: get_RK_a index out of range.\n");
      return aa[ii * ss + jj];
    }

    double get_RK_b(const int &ii) const override
    {
      ASSERT((ii >= 0 && ii < ss),  "Error: get_RK_b index out of range.\n");
      return bb[ii];
    }

    double get_RK_c(const int &ii) const override
    {
      ASSERT((ii >= 0 && ii < ss),"Error: get_RK_c index out of range.\n");
      return cc[ii];
    }

    virtual void print_coefficients() const override
    {
      SYS_T::commPrint("Coefficients:\n");
      SYS_T::commPrint("c: \n");
      for (const auto& ci : cc) SYS_T::commPrint("%e \n", ci);
      SYS_T::commPrint("\nb: \n");
      for (const auto& bi : bb) SYS_T::commPrint("%e \n", bi);

      SYS_T::commPrint("a:\n");
      for (int ii = 0; ii < ss; ++ii) 
      {
        for (const auto& aij : std::vector<double>(aa.begin() + ii * ss, aa.begin() + (ii + 1) * ss))
        {
          SYS_T::commPrint("%e ", aij);
        }
        SYS_T::commPrint("\n");
      }
    }

  private:
    std::array<double, 4> cc{};
    std::array<double, 4> bb{};
    std::array<double, 16> aa{};
};

#endif
