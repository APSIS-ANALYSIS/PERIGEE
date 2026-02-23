#ifndef PLOCASSEM_2X2BLOCK_VMS_HYPERELASTICITY_HPP
#define PLOCASSEM_2X2BLOCK_VMS_HYPERELASTICITY_HPP
// ============================================================================
// PLocAssem_2x2Block_VMS_Hyperelasticity.hpp
//
// This is the local assembly for compressible hyperelasticity using tet4 or hex8
// elements and VMS formulation.
//
// Author: Ju Liu
// Date: Jan 2 2022
// ============================================================================
#include "IPLocAssem_2x2Block.hpp"
#include "MaterialModel_Mixed_Elasticity.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "FEAElementFactory.hpp"
#include "QuadPtsFactory.hpp"

class PLocAssem_2x2Block_VMS_Hyperelasticity : public IPLocAssem_2x2Block
{
  public:
    PLocAssem_2x2Block_VMS_Hyperelasticity(
        const FEType &in_type, const int &in_nqp_v, const int &in_nqp_s,
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        std::unique_ptr<MaterialModel_Mixed_Elasticity> in_matmodel );

    virtual ~PLocAssem_2x2Block_VMS_Hyperelasticity();

    virtual int get_dof_0() const {return 3;}

    virtual int get_dof_1() const {return 1;}

    virtual int get_nLocBas_0() const {return nLocBas;}

    virtual int get_nLocBas_1() const {return nLocBas;}

    virtual int get_snLocBas_0() const {return snLocBas;}

    virtual int get_snLocBas_1() const {return snLocBas;}

    virtual int get_nqpv() const {return nqpv;}  

    virtual int get_nqps() const {return nqps;}

    virtual void Zero_Tangent_Residual();

    virtual void Zero_Residual();

    virtual void Zero_sur_Residual()
    {
      for(int ii=0; ii<sur_size_0; ++ii) sur_Residual0[ii] = 0.0;
    }

    virtual void Assem_Estimate();

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &dot_disp,
        const double * const &dot_velo,
        const double * const &dot_pres,
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_disp,
        const double * const &dot_velo,
        const double * const &dot_pres,
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress );
 
    virtual void Assem_Mass_Residual(
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress );

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

  private:
    const FEType elemType;

    const int nqpv, nqps;

    const std::unique_ptr<FEAElement> elementv, elements;

    const std::unique_ptr<IQuadPts> quadv, quads;

    const double rho0, alpha_f, alpha_m, gamma;

    const int nLocBas, snLocBas, vec_size_0, vec_size_1, sur_size_0;
    
    // useful tensors for the material model
    const std::unique_ptr<MaterialModel_Mixed_Elasticity> matmodel;

    void print_info() const;

    std::array<double, 2> get_tau( const double &dt, const double &Jin, const double &dx ) const;

    Vector_3 get_f(const Vector_3 &pt, const double &tt ) const
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    Vector_3 get_traction( const int &ebc_id,
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out) const;
};

#endif
