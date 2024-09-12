#ifndef SI_ROTATION_INFO_HPP
#define SI_ROTATION_INFO_HPP

#include "Vector_3.hpp"

class SI_rotation_info
{
  public:
    SI_rotation_info(const double &angular, const double &in_thred_time, 
      const Vector_3 &point_xyz, const Vector_3 &angular_direc)
    : target_angular_velo(angular), thred_time(in_thred_time), 
      point_rotated(point_xyz), direction_rotated(angular_direc) {}

    ~SI_rotation_info() = default;

    double get_angular_velo(const double &time) const 
    {
      double angular_velo = target_angular_velo;  

      if (time < thred_time && time >= 0.0)
        angular_velo = 0.5 * target_angular_velo * (1 -  std::cos(MATH_T::PI * time / thred_time));

      return angular_velo;
    }

    double get_target_angular_velo(const double &time) const {return target_angular_velo;}

    Vector_3 get_point_rotated() const {return point_rotated;}

    Vector_3 get_direction_rotated() const {return direction_rotated;}

    double get_thred_time() const {return thred_time;}

    double get_rotated_theta(const double &time) const
    {
      double theta = 0.0;

      if(time < thred_time)
        theta = 0.5 * target_angular_velo * (time - thred_time*std::sin(MATH_T::PI * time / thred_time)/MATH_T::PI);
      else
        theta = 0.5 * target_angular_velo * (time - thred_time*std::sin(MATH_T::PI * time / thred_time)/MATH_T::PI) +  target_angular_velo * (time - thred_time);
      
      return theta;
    }
    
  private:
    // Info of rotation axis
    const double target_angular_velo;

    const double thred_time;

    const Vector_3 point_rotated;

    const Vector_3 direction_rotated;    
};

#endif
