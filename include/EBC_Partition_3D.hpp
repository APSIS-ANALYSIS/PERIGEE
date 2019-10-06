#ifndef EBC_PARTITION_3D_HPP
#define EBC_PARTITION_3D_HPP
// ============================================================================
// EBC_Partition_3D.hpp
//
// Elemental Boundary Condition Partition implementation for three-dimensional
// structural meshes.
//
// Author: Ju Liu
// Date: March 24 2016
// ============================================================================
#include "IEBC_Partition.hpp"
#include "IPart.hpp"
#include "IElemBC.hpp"

class EBC_Partition_3D : public IEBC_Partition
{
  public:
    EBC_Partition_3D( const IPart * const &part,
        const IElemBC * const &ebc_list );

    virtual ~EBC_Partition_3D();

    virtual void write_hdf5(const char * FileName) const;

    virtual void print_info() const;

  private:
    const int cpu_rank;

    std::vector<int> LFront_Elem, LBack_Elem, LLeft_Elem, LRight_Elem, LTop_Elem, LBottom_Elem;

    std::vector<int> Num_LBCElem;
};

#endif
