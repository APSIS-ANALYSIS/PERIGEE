#ifndef EXPLICITRK_HEUNRK2P2S_HPP
#define EXPLICITRK_HEUNRK2P2S_HPP
// ==================================================================
// Explicit Runge–Kutta method: Heun’s method
// Two-stage scheme with second-order accuracy on solution
//
// Ref: J.C. Butcher. Numerical methods for ordinary differential equations, 2016
//
// Date: Apr. 14 2025
// Author: Yujie Sun
// ==================================================================
#include "ITimeMethod_RungeKutta.hpp"

class ExplicitRK_HeunRK2p2s final : public ITimeMethod_RungeKutta 
{
  public:
    ExplicitRK_HeunRK2p2s() : ITimeMethod_RungeKutta(2) 
    {
      // c_i values
      cc[0] = 0.0;
      cc[1] = 1.0; 

      // b_i values
      bb[0] = 0.5; 
      bb[1] = 0.5;

      // a_ij values
      aa[0 * ss + 0] = 0.0;  aa[0 * ss + 1] = 0.0;
      aa[1 * ss + 0] = 1.0;  aa[1 * ss + 1] = 0.0;
    }

    ~ExplicitRK_HeunRK2p2s() override = default;

    double get_RK_a(const int &ii, const int &jj) const override
    {
      ASSERT((ii >= 0 && ii < ss && jj >= 0 && jj < ss), "Error: get_RK_a index out of range.\n");
      return aa[ii * ss + jj];
    }

    double get_RK_b(const int &ii) const override
    {
      ASSERT((ii >= 0 && ii < ss), "Error: get_RK_b index out of range.\n");
      return bb[ii];
    }

    double get_RK_c(const int &ii) const override
    {
      ASSERT((ii >= 0 && ii < ss), "Error: get_RK_c index out of range.\n");
      return cc[ii];
    }

    void print_coefficients() const override
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
    std::array<double, 2> cc{};
    std::array<double, 2> bb{};
    std::array<double, 4> aa{};
};

#endif
