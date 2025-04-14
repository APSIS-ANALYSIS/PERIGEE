#ifndef ITIMEMETHOD_RUNGEKUTTA_HPP
#define ITIMEMETHOD_RUNGEKUTTA_HPP
// ============================================================================
// ITimeMethod_RungeKutta.hpp
// ============================================================================
#include <vector>
#include "Sys_Tools.hpp"

class ITimeMethod_RungeKutta 
{
  public:
    // Constructor, initialize RK coefficients
    ITimeMethod_RungeKutta(int steps, int order, bool implicit)
        : ss(steps), mm(order), isImplicit(implicit), 
          cc(steps, 0.0), bb(steps, 0.0), aa(steps * steps, 0.0)
    {}

    virtual ~ITimeMethod_RungeKutta() = default;

    // Print coefficient (for debugging or viewing)
    virtual void printCoefficients() const 
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

    int get_RK_step() const { return ss; }

    int get_RK_order() const { return mm; }

    double get_RK_a(const int &ii, const int &jj) const 
    {
      ASSERT(!(ii >= 0 && ii < ss && jj >= 0 && jj < ss), "Error: get_RK_a index out of range.\n");
      return aa[ii * ss + jj];
    }

    double get_RK_b(const int &ii) const 
    {
      ASSERT(!(ii >= 0 && ii < ss),  "Error: get_RK_b index out of range.\n");
      return bb[ii];
    }

    double get_RK_c(const int &ii) const 
    {
      ASSERT(!(ii >= 0 && ii < ss),"Error: get_RK_c index out of range.\n");
      return cc[ii];
    }

  protected:
    int ss;  // step
    int mm;  // order
    bool isImplicit;  // Is it an implicit RK?

    std::vector<double> cc;  // ci coefficients
    std::vector<double> bb;  // bi coefficients
    std::vector<double> aa;  // aij coefficients

    virtual void setCoefficients() = 0;
};

#endif
