#ifndef QUADPTS_USERDEFINED_TRIANGLE_HPP
#define QUADPTS_USERDEFINED_TRIANGLE_HPP
// ============================================================================
// QuadPts_UserDefined_Triangle.hpp
// ============================================================================
#include "IQuadPts.hpp"

class QuadPts_UserDefined_Triangle final : public IQuadPts
{
  public:
    QuadPts_UserDefined_Triangle() { reset(); }
    
    ~QuadPts_UserDefined_Triangle() override = default;

    int get_dim() const override {return 3;}

    int get_num_quadPts() const override {return 1;}

    double get_qp(const int &comp) const override
    {return qp[comp];}

    double get_qp(const int &ii, const int &comp) const override
    {return qp[3 * ii + comp];}

    void set_qp(const double &xi, const double &eta) override
    { qp = {{ xi,  eta, 1.0 - xi - eta }}; }

    void reset() override
    {
      constexpr double default_value = 0.333333333333333;
      qp = {{default_value, default_value, default_value}};
    }

    bool check_qp_bound() const override
    {
      constexpr double epsilon = 1.0e-4;

      if( qp[0]>=(0-epsilon) && qp[0]<=(1+epsilon) && 
          qp[1]>=(0-epsilon) && qp[1]<=(1+epsilon) && 
          qp[2]>=(0-epsilon) )
        return true;
      else
        return false;
    }

  private:
    // qp : length 3. Stores the r-s-t coordinates of the
    //      quadrature point.
    //      t = 1 - r - s.
    std::array<double,3> qp;
};

#endif
