#ifndef IMATERIALMODEL_ICH_HPP
#define IMATERIALMODEL_ICH_HPP
// ============================================================================
// IMaterialModel_ich.hpp
// ============================================================================
#include "SymmTensor4_3D.hpp"

class IMaterialModel_ich
{
  public:
    IMaterialModel_ich() = default; 

    virtual ~IMaterialModel_ich() = default;

    virtual void print_info() const = 0;

    virtual std::string get_model_name() const = 0;

    virtual SymmTensor2_3D get_PK_2nd( const Tensor2_3D &F ) const = 0;

    virtual double get_energy( const Tensor2_3D &F ) const = 0;

    virtual SymmTensor4_3D get_PK_Stiffness( const Tensor2_3D &F,
       Tensor2_3D &P_ich, SymmTensor2_3D &S_ich ) const = 0;

    // P_ich := F S_ich
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
       Tensor2_3D &P_ich, SymmTensor2_3D &S_ich ) const
    {
      auto AA = get_PK_Stiffness(F, P_ich, S_ich).full();
      AA.MatMult_1(F);
      AA.MatMult_3(F);
      return AA;
    }

    virtual double get_elastic_mu() const
    {
      SYS_T::commPrint("Warning: IMaterialModel_ich::get_elastic_mu() is not implemented. \n");
      return 0.0;
    }

    virtual Vector_3 get_fibre_dir (const int &dir) const
    {
      SYS_T::commPrint("Warning: IMaterialModel_ich::get_fibre_dir() is not implemented. \n");
      return Vector_3();
    }
};

#endif
