#ifndef RUNGE_KUTTA_BUTCHER_HPP
#define RUNGE_KUTTA_BUTCHERK_HPP

#include "Sys_Tools.hpp"
#include <vector>

class Runge_Kutta_Butcher 
{
  public:
    const int ss;  // step
    const int mm;  // order
    const bool isImplicit;  // Is it an implicit RK?

    std::vector<double> cc;  // ci coefficients
    std::vector<double> bb;  // bi coefficients
    std::vector<std::vector<double>> aa;  // aij coefficients

    // Constructor, initialize RK coefficients
    Runge_Kutta_Butcher(const int &steps, const int &order, const bool &flag);

    ~Runge_Kutta_Butcher() = default;

    // Print coefficient (for debugging or viewing)
    void printCoefficients() const;

    double get_RK_a(const int &ii, const int &jj) const
    {
      if(ii >= ss || jj > ss)
      {
        SYS_T::print_fatal("Error: Runge_Kutta_Butcher::get_RK_a: The index is beyond the steps of RK method. \n");
        return 0.0;
      }
      else
        return aa[ii][jj];
    }

    double get_RK_b(const int &ii) const
    {
      if(ii >= ss)
      {
        SYS_T::print_fatal("Error: Runge_Kutta_Butcher::get_RK_b: The index is beyond the steps of RK method. \n");
        return 0.0;
      }
      else
        return bb[ii];
    }

    double get_RK_c(const int &ii) const
    {
      if(ii >= ss)
      {
        SYS_T::print_fatal("Error: Runge_Kutta_Butcher::get_RK_c: The index is beyond the steps of RK method. \n");
        return 0.0;
      }
      else
        return cc[ii];
    }

  private:
    // Set coefficients based on the given number of steps and order
    void setCoefficients();
};

#endif
