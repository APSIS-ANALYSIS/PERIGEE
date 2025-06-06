#ifndef IPLOCASSEM_HPP
#define IPLOCASSEM_HPP
// ============================================================================
// IPLocAssem.hpp
// Interface for parallel local assembly routine.
//
// This is the pure abstract class for different implementation of
// local assembly routines.
//
// The tangent matrix and the residual vector are stored in the array
// of PetscScalar type: Tangent, Residual.
//
// Author: Ju Liu
// Date: Dec. 3 2013
// ============================================================================
#include "FEAElement.hpp"
#include "ALocal_IEN.hpp"
#include "SymmTensor2_3D.hpp"

class IPLocAssem
{
  public:
    IPLocAssem()
    {
      Tangent  = nullptr;
      Residual = nullptr;
      
      sur_Tangent  = nullptr;
      sur_Residual = nullptr;

      Tangent_ss = nullptr;
      Tangent_sr = nullptr;
      Tangent_rs = nullptr;
      Tangent_rr = nullptr;
      Residual_r = nullptr;
      Residual_s = nullptr;
    }

    virtual ~IPLocAssem(){};

    // ------------------------------------------------------------------------
    // Tangent and Residual of volumetric elements 
    // ------------------------------------------------------------------------
    PetscScalar * Tangent;
    
    PetscScalar * Residual;

    // ------------------------------------------------------------------------
    // Tangent and Residual of surface elements 
    // ------------------------------------------------------------------------
    PetscScalar * sur_Tangent;

    PetscScalar * sur_Residual;
    
    // ------------------------------------------------------------------------
    // Tangent and Residual of sliding-interface
    // ------------------------------------------------------------------------
    PetscScalar * Tangent_ss;
    PetscScalar * Tangent_sr;
    PetscScalar * Tangent_rs;
    PetscScalar * Tangent_rr;

    PetscScalar * Residual_s;
    PetscScalar * Residual_r;
    
    // ------------------------------------------------------------------------
    // ! Get degree of freedom of this problem. In segregated algorithms
    //   this dof returns the fully coupled multiphysics problem's dof.
    // ------------------------------------------------------------------------
    virtual int get_dof() const = 0;

    // ------------------------------------------------------------------------
    // ! Get degree of freedom of the matrix. In segregated algorithms, 
    //   this dof returns the actually implicit solver's dof per node.
    //   In fully coupled fashions, this defaults to the get_dof function.
    // ------------------------------------------------------------------------
    virtual int get_dof_mat() const {return get_dof();}

    // ------------------------------------------------------------------------
    // Return the number of local basis
    // ------------------------------------------------------------------------
    virtual int get_nLocBas() const = 0;

    virtual int get_snLocBas() const = 0;

    // ------------------------------------------------------------------------
    // ! Get the number of ebc functions implemented inside this 
    //   local assembly routine
    // ------------------------------------------------------------------------
    virtual int get_num_ebc_fun() const
    {SYS_T::commPrint("Warning: IPLocAssem::get_num_ebc_fun is not implemented. \n");
      return 0;}

    // ------------------------------------------------------------------------
    // ! Assign all values in Tangent matrix 0.0
    //   Call this function before assembly to zero everything in container
    // ------------------------------------------------------------------------
    virtual void Zero_Tangent_Residual() = 0;

    virtual void Zero_sur_Tangent_Residual()
    {
      SYS_T::print_fatal("Error: Zero_sur_Tangent_Residual is not implemented.\n");
    }

    // ------------------------------------------------------------------------
    // ! Assign all values in Residual vector 0.0
    //   Call this function before assembly to zero everything in container
    // ------------------------------------------------------------------------
    virtual void Zero_Residual() = 0;

    virtual void Zero_sur_Residual()
    {
      SYS_T::print_fatal("Error: Zero_sur_Residual is not implemented. \n");
    }

    // ------------------------------------------------------------------------
    // ! Give nonzero pattern of the sparse matrix 
    // ------------------------------------------------------------------------
    virtual void Assem_Estimate() = 0;

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {SYS_T::commPrint("Warning: this Assem_Residual(...) is not implemented. \n");}

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const double * const &mvelo,
        const double * const &mdisp,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {SYS_T::commPrint("Warning: this Assem_Residual(...) is not implemented. \n");}

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const double * const &mvelo,
        const double * const &mdisp,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual(...) is not implemented. \n");}

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual(...) is not implemented. \n");}

    virtual void Assem_Mass_Residual(
        const double * const &vec_b,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {SYS_T::commPrint("Warning: this Assem_Mass_Residual(...) is not implemented. \n");}

    // ------------------------------------------------------------------------
    // Perform Elemental BC surface integration for elemental BC id ebc_id.
    // Based on ebc_id, the traction forcing function will be called accordingly
    // inside the local assembly routine.
    // ------------------------------------------------------------------------
    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {SYS_T::commPrint("Warning: this Assem_Residual_EBC is not implemented.\n");}

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {SYS_T::commPrint("Warning: this Assem_Residual_EBC is not implemented.\n");}

    virtual void Assem_Residual_BackFlowStab(
        const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {SYS_T::commPrint("Warning: this Assem_Residual_BackFlowStab is not implemented.\n");}

    virtual void Assem_Tangent_Residual_BackFlowStab(
        const double &dt,
        const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual_BackFlowStab is not implemented.\n");}

    virtual void Assem_Residual_EBC_Resistance(
        const int &ebc_id, const double &val,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {SYS_T::commPrint("Warning: this Assem_Residual_EBC_Resistance is not implemented.\n");}

    // ------------------------------------------------------------------------
    // ! Get the model parameter 1
    //   This function is used to pass out the parameters appearing in the weak
    //   form, such as the Reynolds number, Capallarity number, etc.
    //   The definition of this function varies depending on the derived class
    // ------------------------------------------------------------------------
    virtual double get_model_para_1() const
    {
      SYS_T::commPrint("Warning: get_model_para_1() is not implemented. \n");
      return 0.0;
    }

    virtual double get_model_para_2() const
    {
      SYS_T::commPrint("Warning: get_model_para_2() is not implemented. \n");
      return 0.0;
    }

    // ------------------------------------------------------------------------
    // This is a function in local assembly that calculates the flow rate.
    // ------------------------------------------------------------------------
    virtual double get_flowrate( const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {
      SYS_T::commPrint("Warning: get_flowrate() is not implemented. \n");
      return 0.0;
    }

    // ------------------------------------------------------------------------
    // This is a function in local assembly that calculates the pressure
    // integrated over surface: int_{Gamma} p dA
    // as well as the area of the surface: int_{Gamma} 1 dA
    // ------------------------------------------------------------------------
    virtual void get_pressure_area( const double * const &vec, 
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        double &pres, double &area )
    {
      SYS_T::commPrint("Warning: get_pressure_area() is not implemented. \n");
    }

    virtual void Assem_Residual_Weak(
        const double &time, const double &dt,
        const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const int &face_id)
    {SYS_T::commPrint("Warning: this Assem_Residual_Weak is not implemented.\n");}

    virtual void Assem_Tangent_Residual_Weak(
        const double &time, const double &dt,
        const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const int &face_id)
    {SYS_T::commPrint("Warning: this Assem_Tangential_Residual_Weak is not implemented.\n");}

    // for ALE_ns
    virtual void Assem_Residual_Weak(
        const double &time, const double &dt,
        const double * const &sol,
        const double * const &local_mvelo,
        const double * const &local_mdisp,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const int &face_id)
    {SYS_T::commPrint("Warning: this Assem_Residual_Weak_Rotated is not implemented.\n");}

    virtual void Assem_Tangent_Residual_Weak(
        const double &time, const double &dt,
        const double * const &sol,
        const double * const &local_mvelo,
        const double * const &local_mdisp,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const int &face_id)
    {SYS_T::commPrint("Warning: this Assem_Tangential_Residual_Weak_Rotated is not implemented.\n");}

    virtual void Assem_Residual_itf_fixed(
        const int &fixed_qua,
        const double &fixed_qw,
        const double &dt,
        const FEAElement * const &fixed_elementv,
        const FEAElement * const &rotated_elementv,
        const double * const &fixed_local_sol, 
        const double * const &rotated_local_sol)
    {SYS_T::commPrint("Warning: this Assem_Residual_itf_fixed is not implemented.\n");}

    virtual void Assem_Residual_itf_rotated(
        const int &rotated_qua,
        const double &rotated_qw,
        const double &dt,
        const FEAElement * const &rotated_elementv,
        const FEAElement * const &fixed_elementv,
        const double * const &rotated_local_sol,
        const double * const &rotated_local_mvelo,
        const double * const &fixed_local_sol)
    {SYS_T::commPrint("Warning: this Assem_Residual_itf_rotated is not implemented.\n");}

    virtual void Assem_Diag_Tangent_Residual_itf_fixed(
        const int &fixed_qua,
        const double &fixed_qw,
        const double &dt,
        const FEAElement * const &fixed_elementv,
        const FEAElement * const &rotated_elementv,
        const double * const &fixed_local_sol, 
        const double * const &rotated_local_sol)
    {SYS_T::commPrint("Warning: this Assem_Diag_Tangent_Residual_itf_fixed is not implemented.\n");}

    virtual void Assem_Diag_Tangent_Residual_itf_rotated(
        const int &rotated_qua,
        const double &rotated_qw,
        const double &dt,
        const FEAElement * const &rotated_elementv,
        const FEAElement * const &fixed_elementv,
        const double * const &rotated_local_sol,
        const double * const &rotated_local_mvelo,
        const double * const &fixed_local_sol)
    {SYS_T::commPrint("Warning: this Assem_Diag_Tangent_Residual_itf_rotated is not implemented.\n");}

    virtual void Assem_Tangent_itf_MF_fixed(
        const int &fixed_qua,
        const double &fixed_qw,
        const double &dt,
        const FEAElement * const &fixed_elementv,
        const FEAElement * const &rotated_elementv,
        const double * const &fixed_local_sol,
        const double * const &rotated_local_sol)
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual_itf_MF_fixed is not implemented.\n");}

    virtual void Assem_Tangent_itf_MF_rotated(
        const int &rotated_qua,
        const double &rotated_qw,
        const double &dt,
        const FEAElement * const &rotated_elementv,
        const FEAElement * const &fixed_elementv,
        const double * const &rotated_local_sol,
        const double * const &fixed_local_sol,
        const double * const &rotated_local_mvelo)
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual_itf_MF_fixed is not implemented.\n");}
};

#endif
