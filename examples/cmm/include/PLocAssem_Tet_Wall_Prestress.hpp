#ifndef PLOCASSEM_TET_WALL_PRESTRESS_HPP
#define PLOCASSEM_TET_WALL_PRESTRESS_HPP
// ==================================================================
// PLocAssem_Tet_Wall_Prestress.hpp
// 
// Parallel Local Assembly routine for prestress on wall surface.
// ==================================================================
#include "IPLocAssem.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_Tet_Wall_Prestress : public IPLocAssem
{
  public:
    PLocAssem_Tet_Wall_Prestress(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nqp, const double &in_rho, 
        const double &in_vis_mu, const double &in_wall_rho, 
        const double &in_nu, const double &in_kappa,
        const int &elemtype = 501 );

    virtual ~PLocAssem_Tet_Wall_Prestress();

    virtual int get_dof() const {return 4;}

    virtual void Zero_Tangent_Residual()
    {
      SYS_T::print_fatal("Error: Zero_Tangent_Residual is not implemented.\n");
    }

    virtual void Zero_sur_Tangent_Residual()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Residual[ii] = 0.0;
      for(int ii=0; ii<sur_size*sur_size; ++ii) sur_Tangent[ii] = 0.0;
    }

    virtual void Zero_Residual()
    {
      SYS_T::print_fatal("Error: Zero_Residual is not implemented.\n");
    }

    virtual void Zero_sur_Residual()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Residual[ii] = 0.0;
    }

    virtual void Assem_Estimate()
    {
      SYS_T::print_fatal("Error: Assem_Estimate is not implemented.\n");
    }

    virtual void Assem_Tangent_Residual_EBC_Wall(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        const double * const &sol_wall_disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &ele_thickness,
        const double * const &ele_youngsmod,
        const double * const &ele_springconst,
        const double * const &ele_dampingconst,
        const double * const &qua_prestress,
        const IQuadPts * const &quad );

    virtual void get_Wall_CauchyStress(
        const double * const &sol_wall_disp,
        const FEAElement * const &element,
        const double * const &ele_youngsmod,
        std::vector<Matrix_3x3> &stress );

  private:
    // Private data
    const double rho0, vis_mu, alpha_f, alpha_m, gamma;

    // wall properties: density, Poisson ratio, shear correction factor (kappa)
    const double rho_w, nu_w, kappa_w; 

    const int face_nqp; // number of quadrature points for wall

    int snLocBas, sur_size;

    // Private functions
    void print_info() const;

    // Return body force acting on the wall domain
    Vector_3 get_fw( const double &x, const double &y, const double &z,
        const double &t ) const
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }
};

#endif
