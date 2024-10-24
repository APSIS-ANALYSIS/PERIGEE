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
#include "SI_rotation_info.hpp"

class PTime_NS_Solver
{
  public:
    PTime_NS_Solver( const std::string &input_name, 
        const int &input_record_freq, 
        const int &input_renew_tang_freq, 
        const double &input_final_time ) : final_time(input_final_time), 
    sol_record_freq(input_record_freq), 
    renew_tang_freq(input_renew_tang_freq), 
    pb_name(input_name) {}

    ~PTime_NS_Solver() = default;

    void print_info() const;

    void TM_NS_GenAlpha(
        const bool &restart_init_assembly_flag,
        PDNSolution * const &sol_base,
        const PDNSolution * const &init_dot_sol,
        const PDNSolution * const &init_sol,
        const PDNSolution * const &init_mdisp,
        const PDNSolution * const &init_mvelo,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
        const ICVFlowRate * const flr_ptr,
        const APart_Node * const &pNode_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_InflowBC * const &infnbc_part,
        const ALocal_RotatedBC * const &rotnbc_part,
        const ALocal_EBC * const &ebc_part,
        IGenBC * const &gbc,
        const ALocal_WeakBC * const &wbc_part,
        const ALocal_Interface * const &itf_part,
        const SI_rotation_info * const &si_ptr,
        SI_T::SI_solution * const &SI_sol,
        SI_T::SI_quad_point * const &SI_qp,
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
        PNonlinear_NS_Solver * const &nsolver_ptr,
        Mat &shell ) const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    std::string generateNumericSuffix(const int &counter) const 
    { return std::to_string(900000000 + counter); }

    std::string Name_Generator( const int &counter ) const
    { return pb_name + generateNumericSuffix(counter); }
    
    std::string Name_dot_Generator( const int &counter ) const
    { return "dot_" + pb_name + generateNumericSuffix(counter); }

    std::string Name_disp_Generator( const int &counter ) const
    { return "DISP_" + generateNumericSuffix(counter); }

    std::string Name_mvelo_Generator( const int &counter ) const
    { return "MVELO_" + generateNumericSuffix(counter); }
    
    void Write_restart_file(const PDNTimeStep * const &timeinfo,
        const std::string &solname ) const;
    
    //This func may be written into the SI_tools?
    Vector_3 get_radius(const Vector_3 &coor, const SI_rotation_info * const &si_ptr) const 
    { 
      const Vector_3 direction_rotated = si_ptr->get_direction_rotated();
      
      const Vector_3 point_rotated = si_ptr->get_point_rotated();  

      // The vector from the rotation point to the input point
      const Vector_3 point_rotated_to_coor (coor.x() - point_rotated.x(), coor.y() - point_rotated.y(), coor.z() - point_rotated.z());

      const double projectd_length = Vec3::dot_product(point_rotated_to_coor, direction_rotated);
      
      // The projection point of the input point on the rotation axis
      const Vector_3 point_projected (point_rotated.x() +  projectd_length * direction_rotated.x(), point_rotated.y() +  projectd_length * direction_rotated.y(), point_rotated.z() +  projectd_length * direction_rotated.z());
      
      // The vector from the projection point to the input point
      return Vector_3 (coor.x()- point_projected.x(), coor.y()- point_projected.y(), coor.z()- point_projected.z());
    } 

    // Get the current point coordinates for the case of rotation around x/y/z-axis
    Vector_3 get_currPts( const Vector_3 init_pt_xyz,
        const double &tt,
        const SI_rotation_info * const &si_ptr,
        const int &type) const
    {       
      Vector_3 curr_pt_xyz(0, 0, 0);

      const Vector_3 radius_pt = get_radius(init_pt_xyz, si_ptr);

      const double rr = radius_pt.norm2();
      
      double angle = 0.0;

      //rotation axis: case 0: x-axis, case 1: y-axis, case 2: z-axis
      switch(type) 
      {
        case 0:
          angle = MATH_T::get_angle_2d(init_pt_xyz.y(), init_pt_xyz.z());
          angle += si_ptr->get_rotated_theta(tt);
          curr_pt_xyz.x() = init_pt_xyz.x();
          curr_pt_xyz.y() = std::cos(angle) * rr;
          curr_pt_xyz.z() = std::sin(angle) * rr;       
          break;
        case 1: 
          angle = MATH_T::get_angle_2d(init_pt_xyz.z(), init_pt_xyz.x());        
          angle += si_ptr->get_rotated_theta(tt);
          curr_pt_xyz.x() = std::sin(angle) * rr;
          curr_pt_xyz.y() = init_pt_xyz.y();
          curr_pt_xyz.z() = std::cos(angle) * rr;            
          break;            
        case 2: 
          angle = MATH_T::get_angle_2d(init_pt_xyz.x(), init_pt_xyz.y());        
          angle += si_ptr->get_rotated_theta(tt);
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

    // Get the current point coordinates for the case of rotation around any axis (i.e., unit rotation vector (a, b, c))
    // Rodrigues's Formula, theta is the rotation angle
    // [ cos(theta) + a*a(1-cos(theta)),    a*b(1-cos(theta)) - c*sin(theta), b*sin(theta) + a*c(1-cos(theta))  ]
    // [ c*sin(theta) + a*b(1-cos(theta)),  cos(theta) + b*b(1-cos(theta)),   -a*sin(theta) + b*c(1-cos(theta)) ]
    // [ -b*sin(theta) + a*c(1-cos(theta)), a*sin(theta) + b*c(1-cos(theta)), cos(theta) + c*c(1-cos(theta))    ]
    Vector_3 get_currPts(const Vector_3 init_pt_xyz,
        const double &tt,
        const SI_rotation_info * const &si_ptr) const
    {      
      const Vector_3 direction_rotated = si_ptr->get_direction_rotated();
      const Vector_3 point_rotated = si_ptr->get_point_rotated(); 

      const double aa = direction_rotated.x();
      const double bb = direction_rotated.y();       
      const double cc = direction_rotated.z();

      const double theta = si_ptr->get_rotated_theta(tt); 

      const double m00 = std::cos(theta) + aa*aa*(1-std::cos(theta)); 
      const double m01 = aa*bb*(1-std::cos(theta)) - cc*std::sin(theta);
      const double m02 = bb*std::sin(theta) + aa*cc*(1-std::cos(theta));

      const double m10 = cc*std::sin(theta) + aa*bb*(1-std::cos(theta));        
      const double m11 = std::cos(theta) + bb*bb*(1-std::cos(theta)); 
      const double m12 = -aa*std::sin(theta) + bb*cc*(1-std::cos(theta));         

      const double m20 = -bb*std::sin(theta) + aa*cc*(1-std::cos(theta));
      const double m21 = aa*std::sin(theta) + bb*cc*(1-std::cos(theta));
      const double m22 = std::cos(theta) + cc*cc*(1-std::cos(theta));

      const Tensor2_3D mat_rotation (m00, m01, m02, m10, m11, m12, m20, m21, m22);
      
      Vector_3 temp_pt_xyz = init_pt_xyz;    

      temp_pt_xyz -= point_rotated;

      Vector_3 curr_pt_xyz = mat_rotation.VecMult(temp_pt_xyz);

      curr_pt_xyz += point_rotated;

      return curr_pt_xyz;
    }
};

#endif
