#ifndef EXPLICITRK_SSPRK3P3S_HPP
#define EXPLICITRK_SSPRK3P3S_HPP
// ==================================================================
// Explicit Runge–Kutta method: strong-stability preserving method
// Three-stage scheme with third-order accuracy on solution
//
// Ref: S. Gottlieb et al. Strong stability-preserving high-order time 
//      discretization methods. SIAM Review, 2001
// ==================================================================
#include "ITimeMethod_RungeKutta.hpp"

class ExplicitRK_SSPRK3p3s final : public ITimeMethod_RungeKutta 
{
  public:
    ExplicitRK_SSPRK3p3s() : ITimeMethod_RungeKutta(3) 
    {
      // c_i values
      cc[0] = 0.0;
      cc[1] = 1.0;
      cc[2] = 0.5;

      // b_i values
      bb[0] = 1.0 / 6.0;
      bb[1] = 1.0 / 6.0;
      bb[2] = 2.0 / 3.0;

      // a_ij values
      aa[0 * ss + 0] = 0.0;  aa[0 * ss + 1] = 0.0;  aa[0 * ss + 2] = 0.0;
      aa[1 * ss + 0] = 1.0;  aa[1 * ss + 1] = 0.0;  aa[1 * ss + 2] = 0.0;
      aa[2 * ss + 0] = 0.25; aa[2 * ss + 1] = 0.25; aa[2 * ss + 2] = 0.0;
    }

    ~ExplicitRK_SSPRK3p3s() override = default;

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
    std::array<double, 3> cc{};
    std::array<double, 3> bb{};
    std::array<double, 9> aa{};
};

#endif
