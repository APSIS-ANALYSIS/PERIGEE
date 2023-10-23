#ifndef PLOCASSEM_TET_VMS_NS_GENALPHA_HPP
#define PLOCASSEM_TET_VMS_NS_GENALPHA_HPP
// ==================================================================
// PLocAssem_Tet_VMS_NS_GenAlpha.hpp
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

class PLocAssem_Tet_VMS_NS_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_Tet_VMS_NS_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas, const double &in_rho, 
        const double &in_vis_mu, const double &in_beta,
        const double &in_ct = 4.0, const double &in_ctauc = 1.0, 
        const int &elemtype = 501 );

    virtual ~PLocAssem_Tet_VMS_NS_GenAlpha();

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

    virtual void Assem_Residual_Weak(
        const int &weakbc_id,
        const double &time, const double &dt,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const double * const &veleCtrlPts_x,
        const double * const &veleCtrlPts_y,
        const double * const &veleCtrlPts_z,
        const IQuadPts * const &quadv,
        const IQuadPts * const &quads);

    virtual void Assem_Tangential_Residual_Weak(
        const int &weakbc_id,
        const double &time, const double &dt,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const double * const &veleCtrlPts_x,
        const double * const &veleCtrlPts_y,
        const double * const &veleCtrlPts_z,
        const IQuadPts * const &quadv,
        const IQuadPts * const &quads);

  private:
    // Private data
    const double rho0, vis_mu, alpha_f, alpha_m, gamma, beta;

    const int nqp; // number of quadrature points

    double CI; // Constants for stabilization parameters
    
    const double CT; // Constants for stabilization parameters

    const double Ctauc; // Constant scaling factor for tau_C

    int nLocBas, snLocBas, vec_size, sur_size;

    // Private functions
    void print_info() const;

    SymmTensor2_3D get_metric( const std::array<double, 9> &dxi_dx ) const;

    // Return tau_m and tau_c in RB-VMS
    std::array<double, 2> get_tau( const double &dt, 
        const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    // Return tau_bar := (v' G v')^-0.5 x rho0, 
    //        which scales like Time x Density
    // Users can refer to Int. J. Numer. Meth. Fluids 2001; 35: 93–116 for more details
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

    typedef Vector_3 ( PLocAssem_Tet_VMS_NS_GenAlpha::*locassem_tet_vms_ns_funs )( 
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const;

    locassem_tet_vms_ns_funs * flist;

    Vector_3 get_ebc_fun( const int &ebc_id, const Vector_3 &pt, 
        const double &tt, const Vector_3 &n_out ) const
    {
      return ((*this).*(flist[ebc_id]))(pt, tt, n_out);
    }

    Vector_3 get_g_weak(const int &node_idx, const double &tt)
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    // ----------------------------------------------------------------
    // ! get_tau_B : Calculate the coefficient tau_B := [u*]^2 / ||u_tan||,
    //               by solving the non-linear equation of [u+]:
    // g([u+]) = [u+] + 0.1108 * (exp(0.4*[u+]) - 1 - 0.4*[u+] - (0.4*[u+])^2 / 2 - (0.4*[u+])^3 / 6) - [y+]
    //         = 0,
    //               where [u+] := ||u_tan|| / [u*],
    //                     [y+] := y * [u*] / mu = y * ||u_tan|| / (mu * [u+]),
    //               according to Spalding's paper in 1961.
    // Input: \para u_tan : the tangential velocity vector relative to the wall.
    //        \para yy    : the distance from the wall i.e. the 'y' in the formulation.
    //        \para fl_mu : the fluid viscosity i.e the 'mu' in the formulation.
    // ----------------------------------------------------------------
    double get_tau_B(const Vector_3 &u_tan, const double &yy, const double &fl_mu)
    {
      // Use Newton-Raphson method to solve g([u+]) = 0.
      // When [u+] > 0 and [y+] > 0, g([u+]) is monotonically increasing, and there is a unique root.
      const double u_t = u_tan.norm2();

      double u_p0 = 0.0;  // [u+]_i

      double u_p = 1.0;   // [u+]_(i+1)

      do
      {
        u_p0 = u_p;

        const double g_0 = u_p0
        + 0.110803158362334 * (std::exp(0.4*u_p0) - 1.0 - 0.4*u_p0 - 0.08*u_p0*u_p0 - 0.032*u_p0*u_p0*u_p0* (1.0/3.0))
        - yy * u_t * (1.0 / (fl_mu * u_p0));    // g([u+]_i)

        const double g_der_0 = 1
        + 0.110803158362334 * (0.4 * std::exp(0.4*u_p0) - 0.4 - 0.16*u_p0 - 0.032*u_p0*u_p0)
        + yy * u_t * (1.0 / (fl_mu * u_p0 * u_p0)); // dg/d[u+] at [u+]_i

        u_p = u_p0 - g_0 * (1.0 / g_der_0);

      } while (std::abs(u_p - u_p0) > 1.0e-13);
      
      return u_t * (1.0 / (u_p * u_p));   //  tau_B = [u*]^2 / ||u_tan|| = ||u_tan|| / [u+]^2
    }

    // ----------------------------------------------------------------
    // ! build_face_ctrlpt : Given volume element and face id, get face's
    //                     : node coordinates
    // Input: \para ele_type : volume element type
    //        \para face_id  : the face id (Defined in ElemBC_3D::resetSurIEN_outwardnormal)
    //        \para vctrl_.  : the ctrl points of volume element
    // Output: the control points of surface element, corresponding to QuadPts_on_face
    //         and ElemBC_3D::resetSurIEN_outwardnormal
    // ----------------------------------------------------------------
    std::array<std::vector<double>, 3> build_face_ctrlpt( const int &ele_type, const int &face_id,
      const double * const &vctrl_x, const double * const &vctrl_y, const double * const &vctrl_z )
    {
      std::vector<double> fctrl_x {}, fctrl_y {}, fctrl_z {};
      switch (ele_type)
      {
        case 501:
        {
          switch (face_id)
          {
          case 0:
          {
            fctrl_x = std::vector<double> {vctrl_x[3], vctrl_x[1], vctrl_x[2]};
            fctrl_y = std::vector<double> {vctrl_y[3], vctrl_y[1], vctrl_y[2]};
            fctrl_z = std::vector<double> {vctrl_z[3], vctrl_z[1], vctrl_z[2]};
          } break;
          case 1:
          {
            fctrl_x = std::vector<double> {vctrl_x[0], vctrl_x[3], vctrl_x[2]};
            fctrl_y = std::vector<double> {vctrl_y[0], vctrl_y[3], vctrl_y[2]};
            fctrl_z = std::vector<double> {vctrl_z[0], vctrl_z[3], vctrl_z[2]};
          } break;
          case 2:
          {
            fctrl_x = std::vector<double> {vctrl_x[0], vctrl_x[1], vctrl_x[3]};
            fctrl_y = std::vector<double> {vctrl_y[0], vctrl_y[1], vctrl_y[3]};
            fctrl_z = std::vector<double> {vctrl_z[0], vctrl_z[1], vctrl_z[3]};
          } break;
          case 3:
          {
            fctrl_x = std::vector<double> {vctrl_x[0], vctrl_x[2], vctrl_x[1]};
            fctrl_y = std::vector<double> {vctrl_y[0], vctrl_y[2], vctrl_y[1]};
            fctrl_z = std::vector<double> {vctrl_z[0], vctrl_z[2], vctrl_z[1]};
          } break;
          default:
            SYS_T::print_fatal("Error: build_face_ctrlpt, wrong face id of a volume element.\n");
            break;
          }
        } break;
        case 502:
        {
          switch (face_id)
          {
          case 0:
          {
            fctrl_x = std::vector<double> {vctrl_x[3], vctrl_x[1], vctrl_x[2], vctrl_x[5], vctrl_x[9], vctrl_x[8]};
            fctrl_y = std::vector<double> {vctrl_y[3], vctrl_y[1], vctrl_y[2], vctrl_y[5], vctrl_y[9], vctrl_y[8]};
            fctrl_z = std::vector<double> {vctrl_z[3], vctrl_z[1], vctrl_z[2], vctrl_z[5], vctrl_z[9], vctrl_z[8]};
          } break;
          case 1:
          {
            fctrl_x = std::vector<double> {vctrl_x[0], vctrl_x[3], vctrl_x[2], vctrl_x[7], vctrl_x[9], vctrl_x[6]};
            fctrl_y = std::vector<double> {vctrl_y[0], vctrl_y[3], vctrl_y[2], vctrl_y[7], vctrl_y[9], vctrl_y[6]};
            fctrl_z = std::vector<double> {vctrl_z[0], vctrl_z[3], vctrl_z[2], vctrl_z[7], vctrl_z[9], vctrl_z[6]};
          } break;
          case 2:
          {
            fctrl_x = std::vector<double> {vctrl_x[0], vctrl_x[1], vctrl_x[3], vctrl_x[4], vctrl_x[8], vctrl_x[7]};
            fctrl_y = std::vector<double> {vctrl_y[0], vctrl_y[1], vctrl_y[3], vctrl_y[4], vctrl_y[8], vctrl_y[7]};
            fctrl_z = std::vector<double> {vctrl_z[0], vctrl_z[1], vctrl_z[3], vctrl_z[4], vctrl_z[8], vctrl_z[7]};
          } break;
          case 3:
          {
            fctrl_x = std::vector<double> {vctrl_x[0], vctrl_x[2], vctrl_x[1], vctrl_x[6], vctrl_x[5], vctrl_x[4]};
            fctrl_y = std::vector<double> {vctrl_y[0], vctrl_y[2], vctrl_y[1], vctrl_y[6], vctrl_y[5], vctrl_y[4]};
            fctrl_z = std::vector<double> {vctrl_z[0], vctrl_z[2], vctrl_z[1], vctrl_z[6], vctrl_z[5], vctrl_z[4]};
          } break;
          default:
            SYS_T::print_fatal("Error: build_face_ctrlpt, wrong face id of a volume element.\n");
            break;
          }
        } break;
        case 601:
        {
          ; // Unimplement
        } break;
        case 602:
        {
          ; // Unimplement
        } break;
        default:
          SYS_T::print_fatal("Error: build_face_ctrlpt, unknown element type.");
          break;
      }

      return std::array<std::vector<double>, 3> {fctrl_x, fctrl_y, fctrl_z};
    }
};

#endif
