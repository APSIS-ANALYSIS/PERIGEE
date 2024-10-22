#ifndef PGASSEM_2X2BLOCK_FEM_HPP
#define PGASSEM_2X2BLOCK_FEM_HPP
// ==================================================================
// PGAssem_2x2Block_FEM.hpp
//
// Parallel global assembly for block matrices and vectors using the
// AIJ format in PETSc using FEM.
// 
// In addition to the data in the private section, the base class
// IPGAssem_2x2Block contains: Mat K_00, K_01, K_10, and K_11;
// Vec G_0 and G_1.
// 
// Note: In this implementation, we assume the problem is pressure
// velocity coupling in 3D, meaning the dof number for _0 row is 
// 1 (pressure); the dof number for _1 row is 3 (velocity). 
// 
// Also, the periodic boundary condition is not implemented.
//
// Author: Ju Liu
// Date created: Jan 18 2018
// ==================================================================
#include "IPGAssem_2x2Block.hpp"
#include "ALocal_Elem.hpp"

class PGAssem_2x2Block_FEM : public IPGAssem_2x2Block
{
  public:
    PGAssem_2x2Block_FEM(
        IPLocAssem_2x2Block * const &locassem_ptr,
        AGlobal_Mesh_Info const * const &agmi_ptr,
        ALocal_Elem const * const &alelem_ptr,
        ALocal_IEN const * const &aien_ptr,
        APart_Node const * const &pnode_ptr,
        ALocal_NBC const * const &part_nbc,
        ALocal_EBC const * const &part_ebc );

    virtual ~PGAssem_2x2Block_FEM();

    virtual void Get_res_norms( const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        double &norm_0, double &norm_1 ) const;

    virtual void Assem_nonzero_estimate(
        IPLocAssem_2x2Block * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_NBC * const &nbc_part );

    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    virtual void Assem_residual(
        const PDNSolution * const &dot_sol_0,
        const PDNSolution * const &dot_sol_1,
        const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        const double &curr_time,
        const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol_0,
        const PDNSolution * const &dot_sol_1,
        const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        const double &curr_time,
        const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

  private:
    // In addition to the following data, the base class
    // contains Mat K_00/01/10/11 and Vec G_0/1.

    // Presumbly, dof = 7, dof_mat = 4, dof_mat_0 = 1, dof_mat_1 = 3
    const int nlocElem;
    const int dof, dof_mat, dof_mat_0, dof_mat_1;
    const int nLocBas, num_ebc; 
    const int loc_dof_0, loc_dof_1; // dof_mat_0/1 x nLocBas
    int snLocBas;

    PetscInt * row_index_0, * row_index_1;
    PetscInt * srow_index_0, * srow_index_1;

    // array_sol_0/1 stores the local processor's portion of 
    // the solution sol_0/1 including the values on the ghost nodes.
    // The length is nlocalghonode x dof_mat_0/1 
    double * array_sol_0, * array_sol_1;

    // array_dot_0/1 stores the local processor's portion of the
    // solution sol_dot_0/1 including the values on the ghost nodes.
    // The length is nlocalghonode x dof_mat_0/1
    double * array_dot_0, * array_dot_1;

    // local_sol_0/1 stores the element's nodes' solution values.
    // Its length is nLocBas x dof_mat_0/1
    double * local_sol_0, * local_sol_1;

    // local_dot_0/1 stores the element's nodes' solution sol_dot_0/1
    // values. Its length is nLocBas x dof_mat_0/1
    double * local_dot_0, * local_dot_1;

    // local_sur_0/1 stores the local solution values on the 
    // boundary surfaces.
    // Its length is s(urface)nLocBas x dof_mat_0/1
    double * local_sur_0, * local_sur_1;

    // IEN_e: the element's IEN array
    int * IEN_e;

    // LSIEN: local surface element's IEN array
    int * LSIEN;

    // ectrl_x/y/z : the element's nodal control points.
    // The length is nLocBas.
    double * ectrl_x, * ectrl_y, * ectrl_z;

    // sctrl_x/y/z : the surface element's nodal control points.
    // The length is s(urface)nLocBas.
    double * sctrl_x, * sctrl_y, * sctrl_z;

    // Impose essential boundary conditions on field 1,2,3.
    // Field 0 is assumbed to be pressure, which does not have essential bc.
    void EssBC_K11G1( const ALocal_NBC * const &nbc_part );
    void EssBC_G1( const ALocal_NBC * const &nbc_part );

    // GetLocal funciton extract the values from array_xxx_0/1 to local_xxx_0/1.
    void GetLocal_0( const double * const &array, const int * const &IEN,
        const int &val_locbas, double * const &local_array ) const;

    void GetLocal_1( const double * const &array, const int * const &IEN,
        const int &val_locbas, double * const &local_array ) const;

    void NatBC_G1( const double &curr_time, const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const lien_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

};

#endif
