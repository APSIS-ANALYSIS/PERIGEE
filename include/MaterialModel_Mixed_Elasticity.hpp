#ifndef MATERIALMODEL_MIXED_ELASTICITY_HPP
#define MATERIALMODEL_MIXED_ELASTICITY_HPP
// ============================================================================
// MaterialModel_Mixed_Elasticity.hpp
// ============================================================================
#include "IMaterialModel_vol.hpp"
#include "IMaterialModel_ich.hpp"

class MaterialModel_Mixed_Elasticity
{
  public:
    MaterialModel_Mixed_Elasticity( std::unique_ptr<IMaterialModel_vol> in_vmodel,
       std::unique_ptr<IMaterialModel_ich> in_imodel ) 
      : vmodel(std::move(in_vmodel)), imodel(std::move(in_imodel)) {};

    virtual ~MaterialModel_Mixed_Elasticity() = default;

    void print_info() const
    {
      SYS_T::commPrint("\t  Isochoric part: \n");
      imodel->print_info();
      SYS_T::commPrint("\t  Volumetric part: \n");
      vmodel->print_info();
    }

    std::string get_model_name() const
    {return imodel->get_model_name() + vmodel->get_model_name();}

    SymmTensor2_3D get_PK_2nd( const Tensor2_3D &F ) const 
    {return imodel->get_PK_2nd(F);}

    SymmTensor4_3D get_PK_Stiffness( const Tensor2_3D &F, Tensor2_3D &P_iso ) const
    {return imodel->get_PK_Stiffness(F, P_iso);}

    // P_iso := F S_iso
    Tensor2_3D get_PK_1st( const Tensor2_3D &F ) const
    {return imodel->get_PK_1st(F);}
   
    SymmTensor2_3D get_Cauchy_stress( const Tensor2_3D &F ) const
    {return imodel->get_Cauchy_stress(F);}

    Tensor4_3D get_PK_FFStiffness( const Tensor2_3D &F, Tensor2_3D &P_iso ) const
    {return imodel->get_PK_FFStiffness(F, P_iso);}

    double get_elastic_mu() const {return imodel->get_elastic_mu();}

    double get_rho_0() const {return vmodel->get_rho_0();}
    
    double get_rho( const double &pp ) const {return vmodel->get_rho(pp);}

    double get_drho_dp( const double &pp ) const { return vmodel->get_drho_dp(pp);}

    double get_beta( const double &pp ) const {return vmodel->get_beta(pp);}

    double get_dbeta_dp( const double &pp ) const {return vmodel->get_dbeta_dp(pp);}

    double get_elastic_kappa() const {return vmodel->get_elastic_kappa();}

  private:
    const std::unique_ptr<IMaterialModel_vol> vmodel;
    const std::unique_ptr<IMaterialModel_ich> imodel;
};

#endif
