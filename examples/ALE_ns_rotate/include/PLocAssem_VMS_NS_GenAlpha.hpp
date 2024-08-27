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
        const int &elemtype, const double &angular,
        const Vector_3 &point_xyz, const Vector_3 &angular_direc,
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
        const double * const &dot_sol,
        const double * const &sol,
        const double * const &mvelo,
        const double * const &mdisp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        const double * const &mvelo,
        const double * const &mdisp,
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
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual_BackFlowStab(
        const double &dt,
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

    const double angular_velo;

    const Vector_3 point_rotated;

    // Info of rotation axis
    const Vector_3 direction_rotated;    

    // M matrix for tau_m
    //             mm[0], mm[1], mm[2]
    // M = coef *  mm[3], mm[4], mm[5]
    //             mm[6], mm[7], mm[8]
    const double coef;
    const std::array<double, 9> mm; 

    // Private functions
    virtual void print_info() const;

    virtual SymmTensor2_3D get_metric( const std::array<double, 9> &dxi_dx ) const;

    // Return tau_m and tau_c in RB-VMS
    virtual std::array<double, 2> get_tau( const double &dt, 
        const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    // Return tau_bar := (v' G v')^-0.5 x rho0, 
    //        which scales like Time x Density
    // Users can refer to Int. J. Numer. Meth. Fluids 2001; 35: 93–116 
    // for more details
    virtual double get_DC( const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    virtual Vector_3 get_f(const Vector_3 &pt, const double &tt) const
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    virtual Vector_3 get_H1(const Vector_3 &pt, const double &tt, 
        const Vector_3 &n_out ) const
    {
      const double p0 = 0.0;
      return Vector_3( p0*n_out.x(), p0*n_out.y(), p0*n_out.z() );
    }

    typedef Vector_3 ( PLocAssem_VMS_NS_GenAlpha::*locassem_vms_ns_funs )( 
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const;

    locassem_vms_ns_funs * flist;

    virtual Vector_3 get_ebc_fun( const int &ebc_id, const Vector_3 &pt, 
        const double &tt, const Vector_3 &n_out ) const
    {
      return ((*this).*(flist[ebc_id]))(pt, tt, n_out);
    }

    Vector_3 get_Poiseuille_traction1(const Vector_3 &pt, const double &tt, const Vector_3 &n_out) const
    {
      const double delta_P = 1; // pressure
      const double Length = 1; // tube length
      const double fl_mu = 1; // viscosity

      const double z = pt.z(), y = pt.y();

      Tensor2_3D velo_grad ( 0.0, -delta_P * y / (2 * fl_mu * Length), -delta_P * z / (2 * fl_mu * Length),
                             0.0, 0.0, 0.0, 
                             0.0, 0.0, 0.0 );

      // velo_grad *= time_ratio;


      Tensor2_3D pp ( 1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0, 
                      0.0, 0.0, 1.0 );

      Tensor2_3D velo_grad_T = velo_grad;
      velo_grad_T.transpose();

      Tensor2_3D S = - pp + fl_mu * (velo_grad + velo_grad_T);

      return S * n_out;
    }

    Vector_3 get_Poiseuille_traction0(const Vector_3 &pt, const double &tt, const Vector_3 &n_out) const
    {
      const double delta_P = 1; // pressure
      const double Length = 1; // tube length
      const double fl_mu = 1; // viscosity

      const double z = pt.z(), y = pt.y();

      Tensor2_3D velo_grad ( 0.0, -delta_P * y / (2 * fl_mu * Length), -delta_P * z / (2 * fl_mu * Length),
                             0.0, 0.0, 0.0, 
                             0.0, 0.0, 0.0 );

      Tensor2_3D pp ( 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 
                      0.0, 0.0, 0.0 );

      Tensor2_3D velo_grad_T = velo_grad;
      velo_grad_T.transpose();

      Tensor2_3D S = - pp + fl_mu * (velo_grad + velo_grad_T);

      return S * n_out;
    }

    Vector_3 get_cubic_velo_traction(const Vector_3 &pt, const double &tt, const Vector_3 &n_out) const
    {
      const double Q = 0.03 * MATH_T::PI; // 0.03 * pi, v_ave = 3
      const double R = 0.1;

      const double A = R * R * MATH_T::PI; // area
      const double fl_mu = 4.0e-2;

      // ux = 0, uy = 0, uz = v_max * (1 - r^3 / R^3);
      const double x = pt.x(), y = pt.y();
      const double v_max = (5 * Q/ (3 * A));
      const double ra = std::sqrt(x * x + y * y);

      Tensor2_3D velo_grad ( 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 
                             -3 * v_max * x * ra / (R * R * R),
                             -3 * v_max * y * ra / (R * R * R),
                             0.0);

      Tensor2_3D velo_grad_T = velo_grad;
      velo_grad_T.transpose();

      Tensor2_3D S = velo_grad + velo_grad_T;

      S *= fl_mu;
      
      return S * n_out;
    }

    // Get the radius of rotation
    virtual Vector_3 get_radius (const Vector_3 &coor) const
    { 
      // The vector from the rotation point to the input point
      const Vector_3 point_rotated_to_coor (coor.x() - point_rotated.x(), coor.y() - point_rotated.y(), coor.z() - point_rotated.z());

      const double projectd_length = Vec3::dot_product(point_rotated_to_coor, direction_rotated);
      
      // The projection point of the input point on the rotation axis
      const Vector_3 point_projected (point_rotated.x() +  projectd_length * direction_rotated.x(), point_rotated.y() +  projectd_length * direction_rotated.y(), point_rotated.z() +  projectd_length * direction_rotated.z());
      
      // The vector from the projection point to the input point
      return Vector_3 (coor.x()- point_projected.x(), coor.y()- point_projected.y(), coor.z()- point_projected.z());
    } 

    // Get the current point coordinates
    void get_currPts( const double * const &ept_x,
        const double * const &ept_y,
        const double * const &ept_z,
        const double * const &sol,
        const int &len,
        double * const &currPt_x,
        double * const &currPt_y,
        double * const &currPt_z ) const
    {
      for(int ii=0; ii<len; ++ii)
      {
        currPt_x[ii] = ept_x[ii] + sol[3*ii];
        currPt_y[ii] = ept_y[ii] + sol[3*ii+1];
        currPt_z[ii] = ept_z[ii] + sol[3*ii+2];
      }
    }

    // Get the current point coordinates for the case of rotation around x/y/z-axis
    virtual void get_currPts( const double * const &ept_x,
        const double * const &ept_y,
        const double * const &ept_z,
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
        const Vector_3 radius_ept = get_radius(ept_xyz);

        const double rr = radius_ept.norm2();
        
        double angle = 0.0;

        //case 0: x-axis, case 1: y-axis, case 2: z-axis
        switch(type) 
        {
          case 0:
            mag_angular_velo = angular_velo * direction_rotated.x();
            angle = MATH_T::get_angle_2d(ept_xyz(1), ept_xyz(2));        
            angle += mag_angular_velo * tt;
            currPt_x[ii] = ept_x[ii];
            currPt_y[ii] = std::cos(angle) * rr;
            currPt_z[ii] = std::sin(angle) * rr;            
            break;
          case 1: 
            mag_angular_velo = angular_velo * direction_rotated.y();
            angle = MATH_T::get_angle_2d(ept_xyz(2), ept_xyz(0));        
            angle += mag_angular_velo * tt;
            currPt_x[ii] = std::sin(angle) * rr;
            currPt_y[ii] = ept_y[ii];
            currPt_z[ii] = std::cos(angle) * rr;            
            break;            
          case 2: 
            mag_angular_velo = angular_velo * direction_rotated.z();
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
    virtual void get_currPts( const double * const &ept_x,
        const double * const &ept_y,
        const double * const &ept_z,
        const double &tt,
        double * const &currPt_x,
        double * const &currPt_y,
        double * const &currPt_z ) const
    {
      const double aa = direction_rotated.x();
      const double bb = direction_rotated.y();       
      const double cc = direction_rotated.z();

      const double theta = angular_velo * tt; 
        
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

        ept_xyz -= point_rotated;

        Vector_3 cur_xyz = mat_rotation.VecMult(ept_xyz);

        cur_xyz += point_rotated;

        currPt_x[ii] = cur_xyz.x();
        currPt_y[ii] = cur_xyz.y();
        currPt_z[ii] = cur_xyz.z();
      }
    }
};

#endif
