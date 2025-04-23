#ifndef IPGASSEM_HPP
#define IPGASSEM_HPP
// ============================================================================
// IPGAssem.hpp
// Interface for parallel global assembly.
//
// This class provide an interface for different instantiations of global
// assembly of tangent, mass matrix and residual vector, with different BC
// imposed. 
// 
// The data structure for the matrix and the vector are based on the PETSc Mat
// and Vec objects.
//
// Date: July 31 2015 
// ============================================================================
#include "APart_Node.hpp"
#include "APart_Node_FSI.hpp"
#include "ALocal_Elem.hpp"
#include "AGlobal_Mesh_Info.hpp"
#include "IPLocAssem.hpp"
#include "IPLocAssem_2x2Block.hpp"
#include "FEANode.hpp"
#include "PDNSolution.hpp"
#include "ALocal_NBC.hpp"
#include "ALocal_InflowBC.hpp"
#include "ALocal_EBC.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_WeakBC.hpp"
#include "ALocal_Interface.hpp"
#include "Sliding_Interface_Tools.hpp"
#include "IGenBC.hpp"
#include "Tissue_prestress.hpp"

class IPGAssem
{
  public:
    // ------------------------------------------------------------------------
    // Matrix K is a generic matrix object that mainly serves as a SPARSE matrix
    // that we solve with. It can be used as other purposes as well.     
    // ------------------------------------------------------------------------
    Mat K;

    // ------------------------------------------------------------------------
    // Vector G is a generic vector object.
    // ------------------------------------------------------------------------
    Vec G;

    IPGAssem() = default;

    virtual ~IPGAssem() = default;

    // ------------------------------------------------------------------------
    // ! Flag : Fix nonzero structure
    //          Add or insert in new location is ignored. Set after the first
    //          MatAssemblyEnd().
    // ------------------------------------------------------------------------
    void Fix_nonzero_str()
    {MatSetOption(K, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);}

    // ------------------------------------------------------------------------
    // ! Flag : New allocation error
    //          Add or insert in new locations will generate an error message.
    //          Supports AIJ and BAIJ formats.
    // ------------------------------------------------------------------------
    void Fix_nonzero_err_str()
    {MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);}

    // ------------------------------------------------------------------------
    // ! Flag : Ignore new allocation
    //          Add or insert in a new allocation will NOT generate error.
    // ------------------------------------------------------------------------
    void Release_nonzero_err_str()
    {MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);}

    // ------------------------------------------------------------------------
    // ! Flag : Keep nonzero pattern of the matrix K
    // ------------------------------------------------------------------------
    void Keep_nonzero_pattern()
    {MatSetOption(K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);}

    // ------------------------------------------------------------------------
    // ! Clear K and G to be zero
    // ------------------------------------------------------------------------
    void Clear_KG()
    {
      MatZeroEntries(K);
      VecSet(G, 0.0);
    }

    // ------------------------------------------------------------------------
    // ! Clear G to be zero
    // ------------------------------------------------------------------------
    void Clear_G() {VecSet(G, 0.0);}

    // ------------------------------------------------------------------------
    // ! Print the vector G on screen
    // ------------------------------------------------------------------------
    void Print_G() const {VecView(G, PETSC_VIEWER_STDOUT_WORLD);}

    // ------------------------------------------------------------------------
    // ! Assem_nonzero_estimate : Assembly nonzero estimate matrix for K.
    //                            Insert 1.0 to every possible nonzero locations.
    // ------------------------------------------------------------------------
    virtual void Assem_nonzero_estimate()
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented.\n");}

    virtual void Assem_nonzero_estimate( const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented.\n");}

    // ------------------------------------------------------------------------
    // ! Assem_mass_residual : assembly mass matrix and corresponding residual 
    //                         vector for 3D problems WITHOUT pre-existing 
    //                         cached quadrature info.
    // ------------------------------------------------------------------------
    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &mdisp )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented.\n");}

    virtual void Assem_mass_residual( const PDNSolution * const &sol )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented.\n");}

    virtual void Assem_mass_residual(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented.\n");}

    // ------------------------------------------------------------------------
    // ! Assem_residual : assembly residual vector for 3D problem WITHOUT
    //                    pre-existing cached quadrature info.
    // ------------------------------------------------------------------------
    // PGAssem_FSI
    virtual void Assem_Residual(
        const double &curr_time, const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres,
        const PDNSolution * const &dot_velo_np1,
        const PDNSolution * const &velo_np1,
        const PDNSolution * const &disp_np1,
        const IGenBC * const &gbc )
        {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    // PGAssem_Wall_Prestress
    virtual void Assem_Residual(
        const double &curr_time, const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres )
        {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    // Assembly with interface integral
    virtual void Assem_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &mvelo,
        const PDNSolution * const &mdisp,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const IGenBC * const &gbc )
        {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    virtual void Assem_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const double &curr_time,
        const double &dt )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    // ------------------------------------------------------------------------
    // ! Assem_tangent_residual : assembly tangent matrix and residual vector 
    //                            for 3D problem WITHOUT pre-existing cached 
    //                            quadrature info.
    // ------------------------------------------------------------------------
    // PGAssem_FSI
    virtual void Assem_Tangent_Residual(
        const double &curr_time, const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres,
        const PDNSolution * const &dot_velo_np1,
        const PDNSolution * const &velo_np1,
        const PDNSolution * const &disp_np1,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented.\n");}

    virtual void Assem_Tangent_Residual(
        const double &curr_time, const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented.\n");}

    // Assembly with interface integral
    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &mvelo,
        const PDNSolution * const &mdisp,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const IGenBC * const &gbc )
        {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented.\n");}

    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const double &curr_time,
        const double &dt )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented.\n");}

    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented.\n");}

    // --------------------------------------------------------------
    // Assembly boundary integrals
    // --------------------------------------------------------------
    // Assem_surface_flowrate
    // Performs surface integral to calculate the flow rate on the face
    // with ebc_id.
    // Integration is performed for the ebc_part(ebc_id) domain for
    //        v_x n_x + v_y n_y + v_z n_z
    // and then a MPI_Reduce is called to collect the value over multiple
    // CPUs    
    virtual double Assem_surface_flowrate(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const int &ebc_id ) const
    {
      SYS_T::commPrint("Warning: IPGAssem::Assem_surface_flowrate is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id ) const
    {
      SYS_T::commPrint("Warning: IPGAssem::Assem_surface_flowrate is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &vec,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id ) const
    {
      SYS_T::commPrint("Warning: IPGAssem::Assem_surface_flowrate is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &vec,
        const int &ebc_id ) const
    {
      SYS_T::commPrint("Warning: IPGAssem::Assem_surface_flowrate is not implemented. \n");
      return 0.0;
    }

    // Assem_surface_ave_pressure
    // Performs surface integral to calculated the pressure integrated
    // over the surface as well as the surface area. Return the
    // pressure integration divided by the area to give a surface
    // averaged pressure. Before division, the values will be collected
    // from multiple CPUs by MPI_Reduce.
    //  integration p dGamma / integration 1 dGamma
    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &vec,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id ) const
    {
      SYS_T::commPrint("Warning: Assem_surface_ave_pressure is not implemented.\n");
      return 0.0;
    }

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &vec,
        const int &ebc_id ) const
    {
      SYS_T::commPrint("Warning: Assem_surface_ave_pressure is not implemented.\n");
      return 0.0;
    }

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres,
        const int &ebc_id ) const
    {
      SYS_T::commPrint("Warning: Assem_surface_ave_pressure is not implemented.\n");
      return 0.0;
    }

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id ) const
    {
      SYS_T::commPrint("Warning: Assem_surface_ave_pressure is not implemented.\n");
      return 0.0;
    }

    // Update wall prestress at all surface quadrature points
    virtual void Update_Wall_Prestress(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres ) const
    {SYS_T::commPrint("Warning: Update_Wall_Prestress() is not implemented.\n");}

    virtual void write_prestress_hdf5() const
    {SYS_T::commPrint("Warning: write_prestress_hdf5() is not implemented.\n");}

    virtual void Update_SI_situation(
        const PDNSolution * const &sol,
        const PDNSolution * const &mvelo,
        const PDNSolution * const &mdisp )
    {
      SYS_T::print_fatal("Warning: Update_SI_situation() is not implemented. \n");
    }

    virtual void Update_SI_sol(
        const PDNSolution * const &sol )
    {
      SYS_T::print_fatal("Warning: Update_SI_sol() is not implemented. \n");
    }

    virtual const FEANode * Get_fnode()
    {
      SYS_T::print_fatal("Warning: Get_fnode() is not implemented. \n");
      return {};
    }

    virtual const APart_Node * Get_pnode()
    {
      SYS_T::print_fatal("Warning: Get_pnode() is not implemented. \n");
      return {};
    }

};

#endif
