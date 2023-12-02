#ifndef IMATERIALMODEL_HPP
#define IMATERIALMODEL_HPP
// ============================================================================
// IMaterialModel.hpp
//
// Interface for material constitutive relations.
//
// There are three pure virtual functions that need to be implemented
// for each instantiations: 
//   -- print_info : print the basic information of the class.
//   -- get_PK : return the 1st and 2nd Piola-Kirchhoff stresses.
//   -- get_PK_Stiffness : return the 1st and 2nd Piola-Kirchhoff 
//                         stresses as well as the elasticity tensor.
// 
// Additional useful and optional funcitons:
//   -- get_PK_FFStiffness : return the PK stresses and the A tensor
//                           A = partial^2 Psi / partial F partial F.
//   -- get_Cauchy_stress : return the Cauchy stress
//   -- get_strain_energy : return the strain energy
//   -- get_elastic_xxx() : return elastic material parameters
//
// 
// Isochoric-Dialatational split material definition:
// 1. If one uses the traditional isochoric-dialatational split using
//    the traditional way of defining constitutive relations (e.g. 
//    Holzapfel book Sec. 6.4), the PK and and stiffness can be defined
//    straightforwardly.
// 2. If one uses the new mixed formulation, where the pressure p is 
//    treated as an independent variable, all the get_xxx_stress
//    functions and the get_xx_Stiffness return the isochoric part of 
//    the stresses and the elasticity tensor. The isochoric 2nd PK 
//    stress, isochoric Cauchy stress, and the isochoric elasticity
//    tensor are defined in the Holzapfel book, Eqn. (6.90), Eqn. (6.94),
//    Eqn. (6.168), respectively. In addition to that, the dialatational
//    part is replaced by the definition of the density rho = rho(p), its
//    derivative rho' = drho/ dp; the isothermal compressibility 
//    beta = beta(p), and its derivative beta' = dbeta / dp.
//
// Date: Sept. 15 2016
// Author: Ju Liu 
// Contact: liujuy@gmail.com
// ============================================================================
#include "HDF5_Writer.hpp"
#include "HDF5_Reader.hpp"
#include "Tensor4_3D.hpp"

class IMaterialModel
{
  public:
    IMaterialModel() = default;
    
    virtual ~IMaterialModel() = default;

    virtual void print_info() const = 0;

    virtual std::string get_model_name() const
    {
      const std::string output = "unknown";
      SYS_T::commPrint("Warning: IMaterialModel::get_model_name() is not implemented. \n");
      return output;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5") const
    {
      SYS_T::commPrint("Warning: IMaterialModel::write_hdf5() is not implemented. \n");
    }

    virtual void get_PK(const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S) const = 0;

    virtual void get_PK_Stiffness(const Tensor2_3D &F, Tensor2_3D &P, 
        Tensor2_3D &S, Tensor4_3D &CC) const = 0;

    // ------------------------------------------------------------------------
    // Input: F : deformation gradient
    // Output: P : 1st PK
    //         S : 2nd PK
    //         AA : F_iK F_jL C_KILJ
    // ------------------------------------------------------------------------
    virtual void get_PK_FFStiffness( const Tensor2_3D &F, Tensor2_3D &P, 
        Tensor2_3D &S, Tensor4_3D &AA ) const
    {
      get_PK_Stiffness(F, P, S, AA);
      AA.MatMult_1(F);
      AA.MatMult_3(F);
    }

    // ------------------------------------------------------------------------
    // Input: F : deformation gradient
    // Output: sigma : Cauchy stress tensor
    // ------------------------------------------------------------------------
    virtual Tensor2_3D get_Cauchy_stress( const Tensor2_3D &F ) const
    {
      Tensor2_3D P, S;
      get_PK(F, P, S);
      const Tensor2_3D Ft = Ten2::transpose( F );
      Tensor2_3D sigma = P * Ft;
      sigma.scale( (1.0/F.det()) );
      return sigma;
    }

    // ------------------------------------------------------------------------
    // Input: F : deformation gradient
    // Output: sigma : Cauchy stress tensor
    //         aa : J^{-1} F_iI F_jJ F_kK F_lL C_IJKL
    // ------------------------------------------------------------------------
    virtual void get_Cauchy_stiffness( const Tensor2_3D &F, Tensor2_3D &sigma,
       Tensor4_3D &aa ) const
    {
      const double invJ = 1.0 / F.det();
      Tensor2_3D P, S;
      get_PK_Stiffness(F, P, S, aa);
      Tensor2_3D Ft(F); Ft.transpose();
      sigma.MatMult(P, Ft);
      sigma.scale( invJ );

      // Push forward the stiffness tensor
      aa.MatMult_1(F); aa.MatMult_2(F); aa.MatMult_3(F); aa.MatMult_4(F);
      aa.scale( invJ );
    }

    virtual double get_strain_energy( const Tensor2_3D &F ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel::get_strain_energy() is not implemented. \n");
      return 0.0;
    }
    
    // ------------------------------------------------------------------------
    // Interfaces for getting material property parameters
    // ------------------------------------------------------------------------
    virtual double get_elastic_E() const
    {
      SYS_T::commPrint("Warning: IMaterialModel::get_elastic_E() is not implemented. \n");
      return 0.0;
    }

    virtual double get_elastic_rho0() const
    {
      SYS_T::commPrint("Warning: IMaterialModel::get_elastic_rho0() is not implemented. \n");
      return 0.0;
    } 

    virtual double get_elastic_nu() const
    {
      SYS_T::commPrint("Warning: IMaterialModel::get_elastic_nu() is not implemented. \n");
      return 0.0;
    } 

    virtual double get_elastic_lambda() const
    {
      SYS_T::commPrint("Warning: IMaterialModel::get_elastic_lambda() is not implemented. \n");
      return 0.0;
    } 

    virtual double get_elastic_mu() const
    {
      SYS_T::commPrint("Warning: IMaterialModel::get_elastic_mu() is not implemented. \n");
      return 0.0;
    } 
    
    virtual double get_elastic_kappa() const
    {
      SYS_T::commPrint("Warning: IMaterialModel::get_elastic_kappa() is not implemented. \n");
      return 0.0;
    } 

    // ------------------------------------------------------------------------
    // Dialatational properties: density rho, isothermal compressibility beta,
    // and beta's derivative w.r.t. pressure.
    // ------------------------------------------------------------------------
    virtual double get_rho( const double &p ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel::get_rho(p) is not implemented. \n");
      return 0.0;
    } 
    
    virtual double get_drho_dp( const double &p ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel::get_drho_dp(p) is not implemented. \n");
      return 0.0;
    } 

    virtual double get_beta( const double &p ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel::get_beta(p) is not implemented. \n");
      return 0.0;
    }

    virtual double get_dbeta_dp( const double &p ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel::get_dbeta_dp(p) is not implemented. \n");
      return 0.0;
    }

    virtual Vector_3 get_fibre_dir( const int &dir ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel::get_fibre_dir() is not implemented. \n");
      return Vector_3();
    }

    virtual void update_fibre_dir( const Vector_3 &basis_r, const Vector_3 &basis_c, const Vector_3 &basis_l )
    {
      SYS_T::commPrint("Warning: IMaterialModel::update_fibre_dir() is not implemented. \n");
    }
};

#endif
