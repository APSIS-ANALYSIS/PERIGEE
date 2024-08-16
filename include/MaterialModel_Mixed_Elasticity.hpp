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
    IMaterialModel_mixed( std::unique_ptr<IMaterialModel_vol> in_vmodel ) 
      : vmodel(std::move(in_vmodel)) {};

    virtual ~IMaterialModel_mixed() = default;

    virtual void print_info() const = 0;

    virtual std::string get_model_name() const
    {
      SYS_T::commPrint("Warning: IMaterialModel_mixed::get_model_name() is not implemented. \n");
      return std::string("unknown");
    }

    virtual SymmTensor2_3D get_PK_2nd( const Tensor2_3D &F ) const = 0;

    virtual SymmTensor4_3D get_PK_Stiffness( const Tensor2_3D &F,
       Tensor2_3D &P_iso ) const = 0;

    // P_iso := F S_iso
    virtual Tensor2_3D get_PK_1st( const Tensor2_3D &F ) const
    {
      return F * get_PK_2nd(F);
    }
   
    virtual SymmTensor2_3D get_Cauchy_stress( const Tensor2_3D &F ) const
    {
      auto out = get_PK_2nd(F);
      out.push_forward_stress(F);
      return (1.0/F.det()) * out;
    }

    virtual Tensor4_3D get_PK_FFStiffness( const Tensor2_3D &F,
       Tensor2_3D &P_iso ) const
    {
      auto AA = get_PK_Stiffness(F, P_iso).full();
      AA.MatMult_1(F);
      AA.MatMult_3(F);
      return AA;
    }

    virtual double get_elastic_mu() const
    {
      SYS_T::commPrint("Warning: IMaterialModel_mixed::get_elastic_mu() is not implemented. \n");
      return 0.0;
    }

    double get_rho_0() const {return vmodel->get_rho_0();}
    
    double get_rho( const double &pp ) const {return vmodel->get_rho(pp);}

    double get_drho_dp( const double &pp ) const { return vmodel->get_drho_dp(pp);}

    double get_beta( const double &pp ) const {return vmodel->get_beta(pp);}

    double get_dbeta_dp( const double &pp ) const {return vmodel->get_dbeta_dp(pp);}

    double get_elastic_kappa() const {return vmodel->get_elastic_kappa();}

  protected:
    std::unique_ptr<IMaterialModel_vol> vmodel;
};

#endif
