#ifndef PTIME_NS_SOLVER_HPP
#define PTIME_NS_SOLVER_HPP
// ==================================================================
// PTime_NS_Solver.hpp
//
// Parallel time solver for NS.
//
// Author: Ju Liu
// Date: May 23 2017
// ==================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_NS_Solver.hpp"
#include "Sl_rotation_info.hpp"

class PTime_NS_Solver
{
  public:
    PTime_NS_Solver( const std::string &input_name, 
        const int &input_record_freq, const int &input_renew_tang_freq, 
        const double &input_final_time );

    ~PTime_NS_Solver() = default;

    void print_info() const;

    void TM_NS_GenAlpha(
        const bool &restart_init_assembly_flag,
        PDNSolution * const &sol_base,
        const PDNSolution * const &init_dot_sol,
        const PDNSolution * const &init_sol,
        const PDNSolution * const &init_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
        const ICVFlowRate * const flr_ptr,
        const APart_Node * const &pNode_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_InflowBC * const &infnbc_part,
        const ALocal_EBC * const &ebc_part,
        IGenBC * const &gbc,
        const ALocal_WeakBC * const &wbc_part,
        ALocal_Interface * const &itf_part,
        const Sl_rotation_info * const &sl_ptr,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementvs,
        FEAElement * const &elementvs_rotated,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IQuadPts * const &free_quad,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PNonlinear_NS_Solver * const &nsolver_ptr ) const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    std::string Name_Generator( const int &counter ) const;
    
    std::string Name_dot_Generator( const int &counter ) const;

    std::string Name_disp_Generator( const int &counter ) const;
    
    void Write_restart_file(const PDNTimeStep * const &timeinfo,
        const std::string &solname ) const;
    
    //This func may be written into the Sl_tools?
    Vector_3 get_radius(const Vector_3 &coor, const Sl_rotation_info * const &sl_ptr) const 
    { 
      const Vector_3 direction_rotated = sl_ptr->get_direction_rotated();
      
      const Vector_3 point_rotated = sl_ptr->get_point_rotated();  

      // The vector from the rotation point to the input point
      const Vector_3 point_rotated_to_coor (coor.x() - point_rotated.x(), coor.y() - point_rotated.y(), coor.z() - point_rotated.z());

      const double projectd_length = Vec3::dot_product(point_rotated_to_coor, direction_rotated);
      
      // The projection point of the input point on the rotation axis
      const Vector_3 point_projected (point_rotated.x() +  projectd_length * direction_rotated.x(), point_rotated.y() +  projectd_length * direction_rotated.y(), point_rotated.z() +  projectd_length * direction_rotated.z());
      
      // The vector from the projection point to the input point
      return Vector_3 (coor.x()- point_projected.x(), coor.y()- point_projected.y(), coor.z()- point_projected.z());
    } 

    // Get the current point coordinates for the case of rotation around x/y/z-axis
    Vector_3 get_currPts( 
        const Vector_3 init_pt_xyz,
        const double &tt,
        const Sl_rotation_info * const &sl_ptr,
        const int &type) const
    {
      double mag_angular_velo = 0.0; // (rad/s)
      const double angular_velo = sl_ptr->get_angular_velo();
      const Vector_3 direction_rotated = sl_ptr->get_direction_rotated();

      Vector_3 curr_pt_xyz(0, 0, 0);

      const Vector_3 radius_pt = get_radius(init_pt_xyz, sl_ptr);

      const double rr = radius_pt.norm2();
      
      double angle = 0.0;

      //case 0: x-axis, case 1: y-axis, case 2: z-axis
      switch(type) 
      {
        case 0:
          mag_angular_velo = angular_velo * direction_rotated.x();
          angle = MATH_T::get_angle_2d(init_pt_xyz.y(), init_pt_xyz.z());      
          angle += mag_angular_velo * tt;
          curr_pt_xyz.x() = init_pt_xyz.x();
          curr_pt_xyz.y() = std::cos(angle) * rr;
          curr_pt_xyz.z() = std::sin(angle) * rr;       
          break;
        case 1: 
          mag_angular_velo = angular_velo * direction_rotated.y();
          angle = MATH_T::get_angle_2d(init_pt_xyz.z(), init_pt_xyz.x());        
          angle += mag_angular_velo * tt;
          curr_pt_xyz.x() = std::sin(angle) * rr;
          curr_pt_xyz.y() = init_pt_xyz.y();
          curr_pt_xyz.z() = std::cos(angle) * rr;            
          break;            
        case 2: 
          mag_angular_velo = angular_velo * direction_rotated.z();
          angle = MATH_T::get_angle_2d(init_pt_xyz.x(), init_pt_xyz.y());        
          angle += mag_angular_velo * tt;
          curr_pt_xyz.x() = std::cos(angle) * rr;
          curr_pt_xyz.y() = std::sin(angle) * rr;
          curr_pt_xyz.z() = init_pt_xyz.z();            
          break;            
        default:
          SYS_T::print_fatal("Error: PTime_NS_Solver::get_currPts: No such type of rotation axis. \n");
          break;        
      }

      return curr_pt_xyz;
    }
};

#endif
