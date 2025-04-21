#ifndef ITIMEMETHOD_RUNGEKUTTA_HPP
#define ITIMEMETHOD_RUNGEKUTTA_HPP
// ============================================================================
// ITimeMethod_RungeKutta.hpp
// Object: This is the interface for Runge-Kutta method classes.
// It records:
// 1. Number of steps;
//
// Author: Yujie Sun
// Date created: Apr. 14 2025
// ============================================================================
#include <vector>
#include "Sys_Tools.hpp"

class ITimeMethod_RungeKutta 
{
  public:
    // Constructor, initialize the number of steps
    ITimeMethod_RungeKutta(int steps) : ss(steps) {}

    virtual ~ITimeMethod_RungeKutta() = default;

    // Print coefficient (for debugging or viewing)
    virtual void print_coefficients() const 
    {SYS_T::commPrint("Warning: print_coefficients is not implemented. \n");}

    int get_RK_step() const { return ss; }

    // Return the RK coefficients aij
    virtual double get_RK_a(const int &ii, const int &jj) const = 0;

    // Return the RK coefficients bi
    virtual double get_RK_b(const int &ii) const = 0;

    // Return the RK coefficients ci
    virtual double get_RK_c(const int &ii) const = 0;

  protected:
    const int ss;  // num of steps
};

#endif
