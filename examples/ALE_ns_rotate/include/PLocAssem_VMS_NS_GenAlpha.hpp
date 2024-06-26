#ifndef PLOCASSEM_VMS_NS_GENALPHA_HPP
#define PLOCASSEM_VMS_NS_GENALPHA_HPP
// ==================================================================
// PLocAssem_VMS_NS_GenAlpha.hpp
// 
// Parallel Local Assembly routine for VMS and Gen-alpha based NS
// solver.
//
// Author: Ju Liu
// Date: Feb. 10 2020
// ==================================================================
#include "IPLocAssem.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "SymmTensor2_3D.hpp"
#include "Math_Tools.hpp"

class PLocAssem_VMS_NS_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_VMS_NS_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas, const double &in_rho, 
        const double &in_vis_mu, const double &in_beta,
        const int &elemtype,
        const double &in_ct = 4.0, const double &in_ctauc = 1.0 );

    virtual ~PLocAssem_VMS_NS_GenAlpha();

    virtual int get_dof() const {return 4;}

    virtual int get_dof_mat() const {return 4;}

    virtual double get_model_para_1() const {return alpha_f;}

    virtual double get_model_para_2() const {return gamma;}

    virtual void Zero_Tangent_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 0.0;
    }

    virtual void Zero_sur_Tangent_Residual()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Residual[ii] = 0.0;
      for(int ii=0; ii<sur_size*sur_size; ++ii) sur_Tangent[ii] = 0.0;
    }

    virtual void Zero_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
    }

    virtual void Zero_sur_Residual()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Residual[ii] = 0.0;
    }

    virtual void Assem_Estimate()
    {
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 1.0;
    }

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const Vector_3 &angular_velo,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const Vector_3 &angular_velo,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Mass_Residual(
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual double get_flowrate( const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void get_pressure_area( const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad,
        double &pres, double &area );

    virtual void Assem_Residual_EBC_Resistance(
        const int &ebc_id, const double &val,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_BackFlowStab(
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual_BackFlowStab(
        const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

  protected:
    // Private data
    const double rho0, vis_mu, alpha_f, alpha_m, gamma, beta;

    const double CI, CT; // Constants for stabilization parameters
    
    const double Ctauc; // Constant scaling factor for tau_C

    const int nqp; // number of quadrature points
    
    const int nLocBas, snLocBas, vec_size, sur_size;

    // M matrix for tau_m
    //             mm[0], mm[1], mm[2]
    // M = coef *  mm[3], mm[4], mm[5]
    //             mm[6], mm[7], mm[8]
    const double coef;
    const std::array<double, 9> mm; 

    // Private functions
    virtual void print_info() const;

    SymmTensor2_3D get_metric( const std::array<double, 9> &dxi_dx ) const;

    // Return tau_m and tau_c in RB-VMS
    std::array<double, 2> get_tau( const double &dt, 
        const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    // Return tau_bar := (v' G v')^-0.5 x rho0, 
    //        which scales like Time x Density
    // Users can refer to Int. J. Numer. Meth. Fluids 2001; 35: 93–116 
    // for more details
    double get_DC( const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    Vector_3 get_f(const Vector_3 &pt, const double &tt) const
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    Vector_3 get_H1(const Vector_3 &pt, const double &tt, 
        const Vector_3 &n_out ) const
    {
      const double p0 = 0.0;
      return Vector_3( p0*n_out.x(), p0*n_out.y(), p0*n_out.z() );
    }

    typedef Vector_3 ( PLocAssem_VMS_NS_GenAlpha::*locassem_vms_ns_funs )( 
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const;

    locassem_vms_ns_funs * flist;

    Vector_3 get_ebc_fun( const int &ebc_id, const Vector_3 &pt, 
        const double &tt, const Vector_3 &n_out ) const
    {
      return ((*this).*(flist[ebc_id]))(pt, tt, n_out);
    }

    // Get the radius of rotation
    Vector_3 get_radius (const Vector_3 &coor, const Vector_3 &angular_velo) const
    { 
      // Info of rotation axis
      const Vector_3 point_rotated (0.5, 0.0, 0.0);
      Vector_3 direction_rotated (angular_velo);
      
      direction_rotated.normalize();

      // The vector from the rotation point to the input point
      const Vector_3 point_rotated_to_coor (coor.x() - point_rotated.x(), coor.y() - point_rotated.y(), coor.z() - point_rotated.z());

      const double projectd_length = Vec3::dot_product(point_rotated_to_coor, direction_rotated);
      
      // The projection point of the input point on the rotation axis
      const Vector_3 point_projected (point_rotated.x() +  projectd_length * direction_rotated.x(), point_rotated.y() +  projectd_length * direction_rotated.y(), point_rotated.z() +  projectd_length * direction_rotated.z());
      
      // The vector from the projection point to the input point
      return Vector_3 (coor.x()- point_projected.x(), coor.y()- point_projected.y(), coor.z()- point_projected.z());
    } 

    // Get the current point coordinates for the case of rotation around x/y/z-axis
    void get_currPts( const double * const &ept_x,
        const double * const &ept_y,
        const double * const &ept_z,
        const Vector_3 &angular_velo,
        const double &tt,
        double * const &currPt_x,
        double * const &currPt_y,
        double * const &currPt_z,
        const int &type) const
    {
      double mag_angular_velo = 0.0; // (rad/s)

      for(int ii=0; ii<nLocBas; ++ii)
      {
        const Vector_3 ept_xyz (ept_x[ii], ept_y[ii], ept_z[ii]);
        const Vector_3 radius_ept = get_radius(ept_xyz, angular_velo);

        const double rr = radius_ept.norm2();
        
        double angle = 0.0;

         //case 0: x-axis, case 1: y-axis, case 2: z-axis
        switch(type) 
        {
          case 0:
            mag_angular_velo = angular_velo.x();
            angle = MATH_T::get_angle_2d(ept_xyz(1), ept_xyz(2));        
            angle += mag_angular_velo * tt;
            currPt_x[ii] = ept_x[ii];
            currPt_y[ii] = std::cos(angle) * rr;
            currPt_z[ii] = std::sin(angle) * rr;            
            break;
          case 1: 
            mag_angular_velo = angular_velo.y();
            angle = MATH_T::get_angle_2d(ept_xyz(0), ept_xyz(2));        
            angle += mag_angular_velo * tt;
            currPt_x[ii] = std::sin(angle) * rr;
            currPt_y[ii] = ept_y[ii];
            currPt_z[ii] = std::cos(angle) * rr;            
            break;            
          case 2: 
            mag_angular_velo = angular_velo.z();
            angle = MATH_T::get_angle_2d(ept_xyz(0), ept_xyz(1));        
            angle += mag_angular_velo * tt;
            currPt_x[ii] = std::cos(angle) * rr;
            currPt_y[ii] = std::sin(angle) * rr;
            currPt_z[ii] = ept_z[ii];            
            break;            
          default:
            SYS_T::print_fatal("Error: PLocAssem_VMS_NS_GenAlpha::get_currPts: No such type of rotation axis. \n");
            break;        
        }
      }
    }

    // Get the current point coordinates for the case of rotation around any axis (i.e., unit rotation vector (a, b, c))
    // Rodrigues's Formula, theta is the rotation angle
    // [ cos(theta) + a*a(1-cos(theta)),    a*b(1-cos(theta)) - c*sin(theta), b*sin(theta) + a*c(1-cos(theta))  ]
    // [ c*sin(theta) + a*b(1-cos(theta)),  cos(theta) + b*b(1-cos(theta)),   -a*sin(theta) + b*c(1-cos(theta)) ]
    // [ -b*sin(theta) + a*c(1-cos(theta)), a*sin(theta) + b*c(1-cos(theta)), cos(theta) + c*c(1-cos(theta))    ]
    void get_currPts( const double * const &ept_x,
        const double * const &ept_y,
        const double * const &ept_z,
        const Vector_3 &angular_velo,
        const double &tt,
        double * const &currPt_x,
        double * const &currPt_y,
        double * const &currPt_z ) const
    {
      // Info of rotation axis
      const Vector_3 point_rotated (0.5, 0.0, 0.0);

      const double mag_angular_velo = angular_velo.norm2(); // (rad/s)

      // The rotation direction is the direction of angular velocity
      Vector_3 direction_rotated (angular_velo);

      direction_rotated.normalize();

      const double aa= direction_rotated.x();
      const double bb= direction_rotated.y();       
      const double cc= direction_rotated.z();

      const double theta = mag_angular_velo * tt; 
        
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

      for(int ii=0; ii<nLocBas; ++ii)
      {
        Vector_3 ept_xyz (ept_x[ii], ept_y[ii], ept_z[ii]);

        // The vector from the rotation point to the node point
        const Vector_3 point_rotated_to_ept (ept_xyz.x() - point_rotated.x(), ept_xyz.y() - point_rotated.y(), ept_xyz.z() - point_rotated.z());

        const double projectd_length = Vec3::dot_product(point_rotated_to_ept, direction_rotated);
      
        // The projection point of the node point on the rotation axis
        const Vector_3 point_projected (point_rotated.x() +  projectd_length * direction_rotated.x(), point_rotated.y() +  projectd_length * direction_rotated.y(), point_rotated.z() +  projectd_length * direction_rotated.z());
        
        ept_xyz -= point_projected;
        
        Vector_3 cur_xyz = mat_rotation.VecMult(ept_xyz);
        
        cur_xyz += point_projected;

        currPt_x[ii] = cur_xyz.x();
        currPt_y[ii] = cur_xyz.y();
        currPt_z[ii] = cur_xyz.z();
      }
    }
};

#endif
