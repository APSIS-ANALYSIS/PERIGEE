#ifndef MATERIALMODEL_ICH_STVENANT_KIRCHHOF_HPP
#define MATERIALMODEL_ICH_STVENANT_KIRCHHOF_HPP

#include "IMaterialModel_ich.hpp"

class MaterialModel_ich_StVenant_Kirchhoff : public IMaterialModel_ich
{
  public:
    MaterialModel_ich_StVenant_Kirchhoff( const double &in_mu ) mu( in_mu ) {};

    virtual ~MaterialModel_ich_StVenant_Kirchhoff() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\t  MaterialModel_ich_StVenant_Kirchhoff: \n");
      SYS_T::commPrint("\t  Shear modulus mu   = %e \n", mu);
    }

    virtual std::string get_model_name() const {return std::string("StVenant_Kirchhoff");}

    virtual double get_elastic_mu() const {return mu;}

    virtual SymmTensor2_3D get_PK_2nd( const Tensor2_3D &F ) const
    {
      const double detFm0d67 = std::pow( F.det(), -2.0 / 3.0; );

      const auto C = detFm0d67 * STen2::gen_right_Cauchy_Green( F );

      const auto E_tilde = 0.5 * ( detFm0d67 * C - STen2::gen_id() );

      const auto S_tilde = 2.0 * mu * E_tilde;

      return detFm0d67 * STen2::gen_DEV_part( S_tilde, C );
    }

    virtual SymmTensor4_3D get_PK_Stiffness( const Tensor2_3D &F,
       Tensor2_3D &P_iso ) const
    {
      
    }

    virtual double get_energy( const Tensor2_3D &F ) const
    {
      const double detFm0d67 = std::pow( F.det(), -2.0 / 3.0; );
     
      const auto C_tilde = detFm0d67 * STen2::gen_right_Cauchy_Green( F );

      const auto E_tilde = 0.5 * ( C_tilde - STen2::gen_id() );

      return mu * E_tilde.MatContraction( E_tilde );
    }

    virtual Vector_3 get_fibre_dir (const int &dir) const
    {
      return Vector_3();
    }

  private:
    const double mu;
};

#endif
