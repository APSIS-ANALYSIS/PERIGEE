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
#include "IMaterialModel.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_2x2Block_VMS_Hyperelasticity : public IPLocAssem_2x2Block
{
  public:
    PLocAssem_2x2Block_VMS_Hyperelasticity(
        IMaterialModel * const &in_matmodel,
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_snlocbas );

    virtual ~PLocAssem_2x2Block_VMS_Hyperelasticity();

    virtual int get_dof_0() const {return 3;}

    virtual int get_dof_1() const {return 1;}

    virtual int get_nLocBas_0() const {return nLocBas;}

    virtual int get_nLocBas_1() const {return nLocBas;}

    virtual int get_snLocBas_0() const {return snLocBas;}

    virtual int get_snLocBas_1() const {return snLocBas;}

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
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress,
        const IQuadPts * const &quad,
        const std::vector<Vector_3> &eleBasis_r,
        const std::vector<Vector_3> &eleBasis_c,
        const std::vector<Vector_3> &eleBasis_l );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_disp,
        const double * const &dot_velo,
        const double * const &dot_pres,
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress,
        const IQuadPts * const &quad,
        const std::vector<Vector_3> &eleBasis_r,
        const std::vector<Vector_3> &eleBasis_c,
        const std::vector<Vector_3> &eleBasis_l );
 
    virtual void Assem_Mass_Residual(
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress,
        const IQuadPts * const &quad, 
        const std::vector<Vector_3> &eleBasis_r,
        const std::vector<Vector_3> &eleBasis_c,
        const std::vector<Vector_3> &eleBasis_l );

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_Interior_Wall_EBC(
        const double &time,
        const double * const &pres,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    // ------------------------------------------------------------------------
    // This function will calculate the Cauchy stress at every quadrature points
    // within this element. The output stress has length quad -> get_num_quadPts()
    // ------------------------------------------------------------------------
    virtual std::vector<Tensor2_3D> get_Wall_CauchyStress(
        const double * const &disp,
        const double * const &pres,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad ) const;

  private:
    const double rho0, alpha_f, alpha_m, gamma;

    const int nLocBas, snLocBas, vec_size_0, vec_size_1, sur_size_0;
    
    // useful tensors for the material model
    IMaterialModel * const matmodel;

    void print_info() const;

    std::array<double, 2> get_tau( const double &dt, const double &Jin, const double &dx ) const;

    Vector_3 get_f(const Vector_3 &pt, const double &tt ) const
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    // Use pointers to the member functions to facilitate the automatic
    // treatment of ebc surface integration.
    typedef Vector_3 ( PLocAssem_2x2Block_VMS_Hyperelasticity::*locassem_2x2block_vms_comp_ela_fem_funs )( 
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out) const;

    locassem_2x2block_vms_comp_ela_fem_funs * flist;

    Vector_3 get_ebc_fun( const int &ebc_id,
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out) const
    {
      return ((*this).*(flist[ebc_id]))(pt, tt, n_out);
    }
};

#endif
