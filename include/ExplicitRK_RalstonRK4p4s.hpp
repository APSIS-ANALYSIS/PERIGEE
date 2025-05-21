#ifndef EXPLICITRK_RALSTONRK4P4S_HPP
#define EXPLICITRK_RALSTONRK4P4S_HPP
// ==================================================================
// Explicit Runge–Kutta method: Ralston’s method
// Four-stage scheme with fourth-order accuracy on solution
//
// Ref: J.C. Butcher. Numerical methods for ordinary differential equations, 2016
//
// Date: Apr. 14 2025
// Author: Yujie Sun
// ==================================================================
#include "ITimeMethod_RungeKutta.hpp"

class ExplicitRK_RalstonRK4p4s final : public ITimeMethod_RungeKutta 
{
  public:
    ExplicitRK_RalstonRK4p4s() : ITimeMethod_RungeKutta(4) 
    {
      const double beta = std::sqrt(5);        

      // c_i values
      cc[0] = 0.0;
      cc[1] = 2.0 / 5.0;
      cc[2] = (14.0 - 3.0 * beta) / 16.0;
      cc[3] = 1.0;

      // b_i values
      bb[0] = (263.0 + 24.0 * beta) / 1812.0;
      bb[1] = (125.0 - 1000.0 * beta) / 3828.0;
      bb[2] = (3426304.0 + 1661952.0 * beta) / 5924787.0;
      bb[3] = (30.0 - 4.0 * beta) / 123.0;

      // a_ij values
      aa[0 * ss + 0] = 0.0;                                aa[0 * ss + 1] = 0.0;                               aa[0 * ss + 2] = 0.0;                                     aa[0 * ss + 3] = 0.0;
      aa[1 * ss + 0] = 2.0 / 5.0;                          aa[1 * ss + 1] = 0.0;                               aa[1 * ss + 2] = 0.0;                                     aa[1 * ss + 3] = 0.0;
      aa[2 * ss + 0] = (-2889.0 + 1428.0  *beta) / 1024.0; aa[2 * ss + 1] = (3785.0 - 1620.0 * beta) / 1024.0; aa[2 * ss + 2] = 0.0;                                     aa[2 * ss + 3] = 0.0;
      aa[3 * ss + 0] = (-3365.0 + 2094.0 * beta) / 6040.0; aa[3 * ss + 1] = (-975.0 - 3046.0 * beta) / 2552.0; aa[3 * ss + 2] = (467040.0 + 203968.0 * beta) / 240845.0; aa[3 * ss + 3] = 0.0;
    }

    ~ExplicitRK_RalstonRK4p4s() override = default;

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

      SYS_T::commPrint("\na:\n");
      for (int ii = 0; ii < ss; ++ii) 
      {
        for (int jj = 0; jj < ss; ++jj)
          SYS_T::commPrint("%e ", aa[ii * ss + jj]);
      
        SYS_T::commPrint("\n");
      }
    }

  private:
    std::array<double, 4> cc{};
    std::array<double, 4> bb{};
    std::array<double, 16> aa{};
};

#endif
