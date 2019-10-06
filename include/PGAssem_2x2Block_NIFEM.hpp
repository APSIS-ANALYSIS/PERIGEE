#ifndef PGASSEM_2X2BLOCK_NIFEM_HPP
#define PGASSEM_2X2BLOCK_NIFEM_HPP
// ==================================================================
// PGAssem_2x2Block_NIFEM.hpp
//
// Parallel global assembly for block matrices and vectors using the
// AIJ format in PETSc for Non-Isoparametric FEM.
// 
// In the base class, the block data are defined: Mat K_00, K_01,
// K_10, K_11, and Vec G_0 and G_1.
//
// For non-isoparametric finite element implementations, we assume
// that there are two IEN arrays. One is associated with all the
// physics variables; another one is associated with the variable that
// defines the geometry.
// 
// Default setting: the implementation is prepared for mass-linear
// momentum coupling in 3D. For the subscripts, _0 represents quantities
// associated with pressure (scalar); and _1 represents quantities
// associated with velocity (3-vector).
//
// We do NOT consider periodic boundary conditions.
//
// Author: Ju Liu
// Date created: Feb. 7 2018.
// ==================================================================
#include "IPGAssem_2x2Block.hpp"

class PGAssem_2x2Block_NIFEM : public IPGAssem_2x2Block
{
  public:
    // In the constructor, the num_enrichment_nodes gives the number
    // of additional nodes for the physics. The num_sur_enrich_nodes
    // gives the number of additional nodes on surfaces. 
    // For example, low order MINI has one bubble enrichment, 
    // this value is 1 and 0. For P2/P1, there will be P1 enriched 
    // for the physics of pressure, the values are 4 and 3 for linear tet.
    PGAssem_2x2Block_NIFEM(
        const int &in_nlocal_elem,
        const int &in_nlocbas_0,
        const int &in_nlocbas_1,
        const int &in_nlocbas_geo,
        APart_Node const * const &pnode_0_ptr,
        APart_Node const * const &pnode_1_ptr,
        IPLocAssem_2x2Block * const &locassem_ptr,
        ALocal_IEN const * const &aien_0_ptr,
        ALocal_IEN const * const &aien_1_ptr,
        ALocal_IEN const * const &aien_geo_ptr,
        ALocal_NodalBC_2x2Block const * const &part_nbc,
        ALocal_EBC const * const &part_ebc );


    // Destructor
    virtual ~PGAssem_2x2Block_NIFEM();


    // Return norm0 = ||b_0 - K_00 sol_0 - K01 sol1||
    //        norm1 = ||b_1 - k_10 sol_0 - K11 sol1||
    virtual void Get_res_norms( const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        double &norm_0, double &norm_1 ) const;


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
        const ALocal_IEN * const &lien_0_ptr,
        const ALocal_IEN * const &lien_1_ptr,
        const ALocal_IEN * const &lien_geo_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC_2x2Block * const &nbc_part,
        const ALocal_EBC * const &ebc_part );


  private:
    // Base class contains Mat K_00/01/10/11 and Vec G_0/1
    const int nlocElem;
    const int nlocnode_0, nlgn_0; // number of local and local-ghost nodes
    const int nlocnode_1, nlgn_1; // number of local and local-ghost nodes
    const int dof, dof_mat, dof_mat_0, dof_mat_1; // 7,4,1,3
    const int nLocBas_0, nLocBas_1, nLocBas_geo, num_ebc;
    const int loc_dof_0, loc_dof_1; // dof_mat_0/1 x nLocBas_xx
    int snLocBas;

    PetscInt * row_index_0, * row_index_1;
    PetscInt * srow_index_0, * srow_index_1;

    // Stores the local processor's portion of the solution and dot
    // solution, including ghost values. 
    // The length is nlocghonode x dof_mat_0/1
    double * array_sol_0, * array_sol_1;

    double * array_dot_0, * array_dot_1;

    // Stores the element's physical nodes' solution values.
    // Length is nLocBas_phy x dof_mat_0/1
    double * local_sol_0, * local_sol_1;

    double * local_dot_0, * local_dot_1;

    // Stores the solution value on the surface elements
    double * local_sur_0, * local_sur_1;

    // IEN arrays for physics and geometry
    int * IEN_0, * IEN_1, * IEN_geo;

    int * LSIEN;

    // Nodal control points, length is nLocBas_geo
    double * ectrl_x, * ectrl_y, * ectrl_z;

    // Surface control points, length is snLocBas_geo
    double * sctrl_x, * sctrl_y, * sctrl_z;

    // ===== PRIVATE FUNCTIONS =====
    virtual void Assem_nonzero_estimate(
        IPLocAssem_2x2Block * const &lassem_ptr,
        const ALocal_IEN * const &lien_0_ptr,
        const ALocal_IEN * const &lien_1_ptr,
        const ALocal_NodalBC_2x2Block * const &nbc_part );

    void EssBC_K11G1( const ALocal_NodalBC_2x2Block * const &nbc_part );

    void EssBC_G1( const ALocal_NodalBC_2x2Block * const &nbc_part );

    void NatBC_G1( const double &curr_time, const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NodalBC_2x2Block * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    void GetLocal_0( const double * const &array, const int * const &IEN,
        const int &val_locbas, double * const &local_array ) const;

    void GetLocal_1( const double * const &array, const int * const &IEN,
        const int &val_locbas, double * const &local_array ) const;

};


#endif
