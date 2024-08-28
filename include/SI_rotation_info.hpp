#ifndef SI_ROTATION_INFO_HPP
#define SI_ROTATION_INFO_HPP

#include "Vector_3.hpp"

class SI_rotation_info
{
  public:
    SI_rotation_info(const double &angular,const Vector_3 &point_xyz, const Vector_3 &angular_direc)
    : angular_velo(angular), point_rotated(point_xyz), direction_rotated(angular_direc) {}

    ~SI_rotation_info() = default;

    double get_angular_velo() const {return angular_velo;}

    Vector_3 get_point_rotated() const {return point_rotated;}

    Vector_3 get_direction_rotated() const {return direction_rotated;}
    
  private:
    // Info of rotation axis
    const double angular_velo;

    const Vector_3 point_rotated;

    const Vector_3 direction_rotated;    
};

#endif
