#ifndef SL_ROTATION_HPP
#define SL_ROTATION_HPP

#include "Vector_3.hpp"

class Sl_rotation_info
{
  public:
    Sl_rotation_info(const double &angular,const Vector_3 &point_xyz, const Vector_3 &angular_direc)
    : angular_velo(angular), point_rotated(point_xyz), direction_rotated(angular_direc) {}

    ~Sl_rotation_info() = default;

    double get_angular_velo() const {return angular_velo;}

    Vector_3 get_point_rotated() const {return point_rotated;}

    Vector_3 get_direction_rotated() const {return direction_rotated;}

    Vector_3 get_radius(const Vector_3 &coor) const
    { 
      // The vector from the rotation point to the input point
      const Vector_3 point_rotated_to_coor (coor.x() - point_rotated.x(), coor.y() - point_rotated.y(), coor.z() - point_rotated.z());

      const double projectd_length = Vec3::dot_product(point_rotated_to_coor, direction_rotated);
      
      // The projection point of the input point on the rotation axis
      const Vector_3 point_projected (point_rotated.x() +  projectd_length * direction_rotated.x(), point_rotated.y() +  projectd_length * direction_rotated.y(), point_rotated.z() +  projectd_length * direction_rotated.z());
      
      // The vector from the projection point to the input point
      return Vector_3 (coor.x()- point_projected.x(), coor.y()- point_projected.y(), coor.z()- point_projected.z());
    } 

  private:
    // Info of rotation axis
    const double angular_velo;

    const Vector_3 point_rotated;

    const Vector_3 direction_rotated;    
};

#endif
