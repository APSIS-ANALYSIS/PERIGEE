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
    PTime_NS_Solver( std::unique_ptr<PNonlinear_NS_Solver> in_nsolver,
        const std::string &input_name, 
        const int &input_record_freq, 
        const int &input_renew_tang_freq, 
        const double &input_final_time ) : final_time(input_final_time), 
    sol_record_freq(input_record_freq), 
    renew_tang_freq(input_renew_tang_freq), 
    pb_name(input_name), nsolver(std::move(in_nsolver)) {}

    ~PTime_NS_Solver() = default;

    void print_info() const;

    void print_lsolver_info() const {nsolver -> print_lsolver_info();}

    // ------------------------------------------------------------------------
    // Generate a file name for inlet/outlet face as prefix_xxx_data.txt
    // ------------------------------------------------------------------------
    std::string gen_flowfile_name(const std::string &prefix, const int &id) const
    {
      std::ostringstream ss;
      ss << prefix;

      if(id < 10) ss << "00";
      else if(id < 100) ss << "0";

      ss << id << "_data.txt";
      return ss.str();
    }
    
    void TM_NS_GenAlpha(
        const bool &restart_init_assembly_flag,
        std::unique_ptr<PDNSolution> init_dot_sol,
        std::unique_ptr<PDNSolution> init_sol,
        std::unique_ptr<PDNSolution> init_mdisp,
        std::unique_ptr<PDNSolution> init_mvelo,
        std::unique_ptr<PDNTimeStep> time_info,
        std::unique_ptr<ALocal_InflowBC> infnbc_part,
        std::unique_ptr<ALocal_RotatedBC> rotnbc_part,
        std::unique_ptr<IGenBC> gbc,
        std::unique_ptr<SI_rotation_info> rot_info,
        std::unique_ptr<IPGAssem> gassem_ptr,
        Mat &shell ) const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    const std::unique_ptr<PNonlinear_NS_Solver> nsolver;

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
    Vector_3 get_radius(const Vector_3 &coor, const SI_rotation_info * const &rot_info) const 
    { 
      const Vector_3 direction_rotated = rot_info->get_direction_rotated();
      
      const Vector_3 point_rotated = rot_info->get_point_rotated();  

      // The vector from the rotation point to the input point
      const Vector_3 point_rotated_to_coor (coor.x() - point_rotated.x(), coor.y() - point_rotated.y(), coor.z() - point_rotated.z());

      const double projectd_length = Vec3::dot_product(point_rotated_to_coor, direction_rotated);
      
      // The projection point of the input point on the rotation axis
      const Vector_3 point_projected (point_rotated.x() +  projectd_length * direction_rotated.x(), point_rotated.y() +  projectd_length * direction_rotated.y(), point_rotated.z() +  projectd_length * direction_rotated.z());
      
      // The vector from the projection point to the input point
      return Vector_3 (coor.x()- point_projected.x(), coor.y()- point_projected.y(), coor.z()- point_projected.z());
    } 

    // Get the current point coordinates for the case of rotation around any axis (i.e., unit rotation vector (a, b, c))
    // Rodrigues's Formula, theta is the rotation angle
    // [ cos(theta) + a*a(1-cos(theta)),    a*b(1-cos(theta)) - c*sin(theta), b*sin(theta) + a*c(1-cos(theta))  ]
    // [ c*sin(theta) + a*b(1-cos(theta)),  cos(theta) + b*b(1-cos(theta)),   -a*sin(theta) + b*c(1-cos(theta)) ]
    // [ -b*sin(theta) + a*c(1-cos(theta)), a*sin(theta) + b*c(1-cos(theta)), cos(theta) + c*c(1-cos(theta))    ]
    Vector_3 get_currPts(const Vector_3 init_pt_xyz,
        const double &tt,
        const SI_rotation_info * const &rot_info) const
    {      
      const Vector_3 direction_rotated = rot_info->get_direction_rotated();
      const Vector_3 point_rotated = rot_info->get_point_rotated(); 

      const double aa = direction_rotated.x();
      const double bb = direction_rotated.y();       
      const double cc = direction_rotated.z();

      const double theta = rot_info->get_rotated_theta(tt); 

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
