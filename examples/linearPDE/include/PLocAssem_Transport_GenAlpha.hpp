#ifndef PLOCASSEM_TRANSPORT_GENALPHA_HPP
#define PLOCASSEM_TRANSPORT_GENALPHA_HPP
// ============================================================================
// PLocAssem_Transport_GenAlpha.hpp
//
// Parallel Local Assembly routine for Gen-alpha based transport solver。
//
// Date: Jan 21 2022
// ============================================================================
#include "IPLocAssem.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "FEAElementFactory.hpp"
#include "QuadPtsFactory.hpp"

class PLocAssem_Transport_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_Transport_GenAlpha( 
        const FEType &in_type, const int &in_nqp_v, const int &in_nqp_s,
        const double &in_rho, const double &in_cap, const double &in_kappa,
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_num_ebc_fun );

    virtual ~PLocAssem_Transport_GenAlpha();

    virtual int get_dof() const {return 1;}

    virtual int get_dof_mat() const {return 1;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual int get_snLocBas() const {return snLocBas;}

    virtual void Zero_Tangent_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 0.0;
    }

    virtual void Zero_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
    }

    virtual void Zero_sur_Tangent_Residual()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Residual[ii] = 0.0;
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
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    virtual void Assem_Mass_Residual(
        const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

  private:
    // Private data
    const FEType elemType;
    const int nqpv, nqps;
    const std::unique_ptr<FEAElement> elementv, elements;
    const std::unique_ptr<IQuadPts> quadv, quads;

    const double rho, cap, kappa;
    const double alpha_f, alpha_m, gamma;

    const int num_ebc_fun;

    const int nLocBas, snLocBas;
    const int vec_size, sur_size;

    void print_info() const;

    double get_f( const Vector_3 &pt, const double &tt ) const
    {
      //const double pi = MATH_T::PI;

      const double t3 = tt*tt*tt;
      const double t4 = t3 * tt;
      const double x = pt.x();
      const double y = pt.y();
      const double z = pt.z();

      return 4*cap*rho*t3*x*y*z*(x - 1)*(y - 1)*(z - 1) - 2*kappa*t4*x*z*(x - 1)*(z - 1) - 2*kappa*t4*y*z*(y - 1)*(z - 1) - 2*kappa*t4*x*y*(x - 1)*(y - 1); 
    }

    typedef double ( PLocAssem_Transport_GenAlpha::*locassem_transport_funs )( const Vector_3 &pt, const double &t ) const;

    locassem_transport_funs * flist;

    double get_ebc_fun( const int &ebc_id, const Vector_3 &pt, const double &tt ) const
    {
      return ((*this).*(flist[ebc_id]))(pt, tt);
    }

    double get_g_0( const Vector_3 &pt, const double &time ) const
    {
      return 0.0;
    }

    double get_g_1( const Vector_3 &pt, const double &time ) const
    {
      return 0.0;
    }
};

#endif
