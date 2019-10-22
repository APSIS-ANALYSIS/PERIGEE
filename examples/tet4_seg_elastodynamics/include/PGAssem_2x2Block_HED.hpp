#ifndef PGASSEM_2X2BLOCK_HED_HPP
#define PGASSEM_2X2BLOCK_HED_HPP
// ==================================================================
// PGAssem_2x2Block_HED.hpp
//
// Parallel global assembly for block matrices and vectors for
// Hyper-Elasto-Dynamics (HED).
//
// Inherit Mat K_00, K_01, K10, K11 & Vec G_0 & G_1 from the base
// class IPGAssem_2x2Block.
// 
// Author: Ju Liu
// Date created: Feb. 20 2018
// ==================================================================
#include "IPGAssem_2x2Block.hpp"

class PGAssem_2x2Block_HED : public IPGAssem_2x2Block
{
  public:
    PGAssem_2x2Block_HED(
        IPLocAssem_2x2Block * const &locassem_ptr,
        ALocal_IEN const * const &aien_ptr,
        APart_Node const * const &pnode_ptr,
        ALocal_NodalBC const * const &part_nbc,
        ALocal_EBC const * const &part_ebc );

    virtual ~PGAssem_2x2Block_HED();

    virtual void Get_res_norms( const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        double &norm_0, double &norm_1 ) const;

    virtual void Assem_nonzero_estimate(
        IPLocAssem_2x2Block * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_NodalBC * const &nbc_part );

    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_d,
        const PDNSolution * const &sol_p,
        const PDNSolution * const &sol_v,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    virtual void Assem_residual(
        const PDNSolution * const &dot_d,
        const PDNSolution * const &dot_p,
        const PDNSolution * const &dot_v,
        const PDNSolution * const &sol_d,
        const PDNSolution * const &sol_p,
        const PDNSolution * const &sol_v,
        const double &curr_time,
        const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_d,
        const PDNSolution * const &dot_p,
        const PDNSolution * const &dot_v,
        const PDNSolution * const &sol_d,
        const PDNSolution * const &sol_p,
        const PDNSolution * const &sol_v,
        const double &curr_time,
        const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

  private:
    // ** DATA SECTION **
    const int nlocElem;
    const int dof, dof_mat, dof_mat_0, dof_mat_1;
    const int nLocBas, num_ebc;
    const int loc_dof_0, loc_dof_1;
    int snLocBas;

    PetscInt * row_index_0, * row_index_1;
    PetscInt * srow_index_0, * srow_index_1;

    // Local solution vector copy including the value on ghost nodes
    // length is nlocalghonode x dof_mat_0/1 0 for pres, 1 for disp velo
    double * array_sol_d, * array_sol_p, * array_sol_v;
    double * array_dot_d, * array_dot_p, * array_dot_v;

    // Element's copy for solution values.
    // length is nLocBas x dof_mat_0/1
    double * local_sol_d, * local_sol_p, * local_sol_v;
    double * local_dot_d, * local_dot_p, * local_dot_v;

    // Surface element's copy for solution values.
    // length is snLocBas x dof_mat_0/1
    double * local_sur_d, * local_sur_p, * local_sur_v;

    // Element's IEN array
    int * IEN_e;

    // Surface element's IEN array
    int * LSIEN;

    // Element's geometrical control points
    // Length is nLocBas
    double * ectrl_x, * ectrl_y, * ectrl_z;

    // Surface element's geometrical control points
    double * sctrl_x, * sctrl_y, * sctrl_z;


    // ** FUNCTION SECTION **
    void EssBC_K11G1( const ALocal_NodalBC * const &nbc_part );
   

    void EssBC_G1( const ALocal_NodalBC * const &nbc_part );


    void NatBC_G1( const double &curr_time, const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const lien_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );


    void GetLocal_0( const double * const &array, const int * const &IEN,
        const int &val_locbas, double * const &local_array ) const
    {
      for(int ii=0; ii<val_locbas; ++ii) local_array[ii] = array[IEN[ii]];
    }


    void GetLocal_1( const double * const &array, const int * const &IEN,
        const int &val_locbas, double * const &local_array ) const
    {
      for(int ii=0; ii<val_locbas; ++ii)
      {
        local_array[3*ii  ] = array[ 3*IEN[ii] ];
        local_array[3*ii+1] = array[ 3*IEN[ii]+1 ];
        local_array[3*ii+2] = array[ 3*IEN[ii]+2 ];
      }
    }

};


#endif
