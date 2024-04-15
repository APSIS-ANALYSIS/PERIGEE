#ifndef PGASSEM_STRESS_RECOVERY_HPP
#define PGASSEM_STRESS_RECOVERY_HPP

#include "IPGAssem.hpp"
#include "PETSc_Tools.hpp"

class PGAssem_Stress_Recovery : public IPGAssem
{
  public:
    // Constructor
    PGAssem_Stress_Recovery(
        IPLocAssem * const &locassem_ptr,
        const IAGlobal_Mesh_Info * const &agmi_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_ptr,
        const APart_Node * const &pnode_ptr,
        const int &in_nz_estimate = 60);
    
    // Destructor
    virtual ~PGAssem_Stress_Recovery();

    // Nonzero pattern estimate
    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &pnode_ptr );
    
    // Assembly the residual vector
    virtual void Assem_residual(
        const PDNSolution * const &isol,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        const IQuadPts * const &quad_v,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const APart_Node * const &pnode_ptr,
        const Tissue_property * const &tp_ptr );
    
    // Assembly the residual and mass matrix
    virtual void Assem_mass_residual(
        const PDNSolution * const &isol,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        const IQuadPts * const &quad_v,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const APart_Node * const &pnode_ptr,
        const Tissue_property * const &tp_ptr );
    
    private:
      // Private data
      const int nLocBas, nlgn;

      void GetLocal( const double * const &array, const int * const &IEN,
        double * const &local_array ) const
      {
        for(int ii=0; ii<nLocBas; ++ii)
        {
          const int offset1 = ii * 3;
          const int offset2 = IEN[ii] * 3;
          for(int jj=0; jj<3; ++jj)
            local_array[offset1 + jj] = array[offset2 + jj];
          }
      }
};

#endif