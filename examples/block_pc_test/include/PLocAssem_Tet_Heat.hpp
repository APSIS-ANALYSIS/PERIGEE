#ifndef PLOCASSEM_TET_HEAT_HPP
#define PLOCASSEM_TET_HEAT_HPP
// ==================================================================
// PLocAssem_Tet_Heat.hpp
//
// Parallel local assembly routine for heat equation for testing
// laplace operator.
//
// Date: May 12 2020
// ==================================================================
#include "IPLocAssem.hpp"

class PLocAssem_Tet_Heat : public IPLocAssem
{
  public:
    PLocAssem_Tet_Heat( const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas, const int &elemType = 501 );

    virtual ~PLocAssem_Tet_Heat();

    virtual int get_dof() const {return 1;}

    virtual int get_dof_mat() const {return 1;}

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
        const IQuadPts * const &quad )
    {
      SYS_T::print_fatal("Error: not implemented.\n");
    }

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
        const IQuadPts * const &quad )
    {
      SYS_T::print_fatal("Error: not implemented.\n");
    }

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_EBC_Resistance(
        const int &ebc_id, const double &val,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

  private:
    const int nqp;

    int nLocBas, snLocBas, vec_size, sur_size;

    std::vector<double> R, dR_dx, dR_dy, dR_dz;

    void print_info() const;

    void get_f(const double &x, const double &y, const double &z,
        const double &t, double &f ) const
    {
      f = 4.0;
    }

    void get_H1(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &g ) const
    {
      g = 0.0;
    }

    typedef void ( PLocAssem_Tet_Heat::*locassem_tet_heat_funs )( 
        const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &g ) const;

    locassem_tet_heat_funs * flist;

    void get_ebc_fun( const int &ebc_id,
        const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &g ) const
    {
      return ((*this).*(flist[ebc_id]))(x,y,z,t,nx,ny,nz,g);
    }
};

#endif
