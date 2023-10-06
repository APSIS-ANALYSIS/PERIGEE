#ifndef IPLOCASSEM_LINEAR_HPP
#define IPLOCASSEM_LINEAR_HPP
// ============================================================================
// IPLocAssem_Linear.hpp
// Interface for parallel local assembly routine.
//
// This is the pure abstract class for different implementation of
// local assembly routines.
//
// The mass matrix, stiffness matrix, and load vector are stored in the array
// of PetscScalar type: Mass, Stiffness, and Load.
//
// Author: Xinhai Yue
// Date: Oct. 6 2033
// ============================================================================
#include "FEAElement.hpp"
#include "ALocal_IEN.hpp"

class IPLocAssem_Linear
{
  public:
    IPLocAssem_Linear()
    {
      Mass = nullptr;
      Stiffness = nullptr;
      Load = nullptr;
      
      sur_Mass = nullptr;
      sur_Load = nullptr;
    }
    
    virtual ~IPLocAssem_Linear(){};

    // ------------------------------------------------------------------------
    // Mass, Stiffness, Load of volumetric elements 
    // ------------------------------------------------------------------------
    // -------------------------------------------------------------- 
    PetscScalar * Mass;
    PetscScalar * Stiffness
    PetscScalar * Load;

    // ------------------------------------------------------------------------
    // -------------------------------------------------------------- 
    // Mass and Load of surface elements 
    // ------------------------------------------------------------------------
    // -------------------------------------------------------------- 
    PetscScalar * Mass;
    PetscScalar * Load;

    // -------------------------------------------------------------- 
    // ------------------------------------------------------------------------
    // ! Get degree of freedom of this problem. In segregated algorithms
    //   this dof returns the fully coupled multiphysics problem's dof.
    // ------------------------------------------------------------------------
    // -------------------------------------------------------------- 
    virtual int get_dof() const = 0;

    // -------------------------------------------------------------- 
    // ------------------------------------------------------------------------
    // ! Get degree of freedom of the matrix. In segregated algorithms, 
    //   this dof returns the actually implicit solver's dof per node.
    //   In fully coupled fashions, this defaults to the get_dof function.
    // ------------------------------------------------------------------------
    // -------------------------------------------------------------- 
    virtual int get_dof_mat() const {return get_dof();}

    // --------------------------------------------------------------
    // ------------------------------------------------------------------------
    // Return the number of local basis
    // ------------------------------------------------------------------------
    // --------------------------------------------------------------
    virtual int get_nLocBas() const
    {
      SYS_T::commPrint("Warning: IPLocAssem_Linear::get_nLocBas is not implemented. \n");
      return -1;
    }

    virtual int get_snLocBas() const
    {
      SYS_T::commPrint("Warning: IPLocAssem_Linear::get_snLocBas is not implemented. \n");
      return -1;
    }

    // -------------------------------------------------------------- 
    // ------------------------------------------------------------------------
    // ! Get the number of ebc functions implemented inside this 
    //   local assembly routine
    // ------------------------------------------------------------------------
    // -------------------------------------------------------------- 
    virtual int get_num_ebc_fun() const
    {SYS_T::commPrint("Warning: IPLocAssem_Linear::get_num_ebc_fun is not implemented. \n");
      return 0;}

    // -------------------------------------------------------------- 
    // ------------------------------------------------------------------------
    // ! Assign all values in Mass and Stiffness matrices 0.0
    //   Call this function before assembly to zero everything in container
    // ------------------------------------------------------------------------
    // -------------------------------------------------------------- 
    virtual void Zero_Mass_Stiffness_Load() = 0;

    virtual void Zero_sur_Mass_Load()
    {
      SYS_T::print_fatal("Error: Zero_sur_Mass_Load is not implemented.\n");
    }

    // -------------------------------------------------------------- 
    // ------------------------------------------------------------------------
    // ! Assign all values in Load vector 0.0
    //   Call this function before assembly to zero everything in container
    // ------------------------------------------------------------------------
    // -------------------------------------------------------------- 
    virtual void Zero_Load() = 0;

    virtual void Zero_sur_Load()
    {
      SYS_T::print_fatal("Error: Zero_sur_Load is not implemented. \n");
    }

    // -------------------------------------------------------------- 
    // ------------------------------------------------------------------------
    // ! Give nonzero pattern of the sparse matrix 
    // ------------------------------------------------------------------------
    // -------------------------------------------------------------- 
    virtual void Assem_Estimate() = 0;

    // ------------------------------------------------------------------------
    // \para element: the container for classical element routine. It only
    //                requires the x-y-z coordinates for the nodes and the 
    //                volumetric quadrature routine to generate the basis
    //                functions
    // ------------------------------------------------------------------------
    virtual void Assem_Load(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Load(...) is not implemented. \n");}

    // ------------------------------------------------------------------------
    // \para element: the container for classical element routine. It only
    //                requires the x-y-z coordinates for the nodes and the 
    //                volumetric quadrature routine to generate the basis
    //                functions
    // ------------------------------------------------------------------------
    virtual void Assem_Stiffness(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Stiffness(...) is not implemented. \n");}

    // ------------------------------------------------------------------------
    // \para element: the container for classical element routine. It only
    //                requires the x-y-z coordinates for the nodes and the 
    //                volumetric quadrature routine to generate the basis
    //                functions
    // ------------------------------------------------------------------------
    virtual void Assem_Mass(
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Mass(...) is not implemented. \n");}

    // ------------------------------------------------------------------------
    // Perform Elemental BC surface integration for elemental BC id ebc_id.
    // Based on ebc_id, the traction forcing function will be called accordingly
    // inside the local assembly routine.
    // \para element: the container for the element, only requires the geometry
    //                information for the control points, and the quadrature
    //                info to generate basis function info.
    // ------------------------------------------------------------------------
    virtual void Assem_Load_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Load_EBC is not implemented.\n");}

    virtual void Assem_Load_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Load_EBC is not implemented.\n");}

    // ------------------------------------------------------------------------
    // Perform Elemental BC surface integration for elemental BC id ebc_id.
    // Based on ebc_id, the traction forcing function will be called accordingly
    // inside the local assembly routine.
    // \para in_x, in_y, in_z : interior point to the element
    // \para element: the container for the element, only requires the geometry
    //                information for the control points, and the quadrature
    //                info to generate basis function info.
    // ------------------------------------------------------------------------
    virtual void Assem_Load_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double &in_x, const double &in_y, const double &in_z,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Load_EBC is not implemented.\n");}
};

#endif
