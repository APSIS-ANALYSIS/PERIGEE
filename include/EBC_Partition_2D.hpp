#ifndef EBC_PARTITION_2D_HPP
#define EBC_PARTITION_2D_HPP
// ============================================================================
// EBC_Partition_2D.hpp
//
// Elemental Boundary Condition Partition implementation for two-dimensional
// structural meshes.
//
// Date: Aug. 19 2015
// ============================================================================
#include "IEBC_Partition.hpp"
#include "IPart.hpp"
#include "IElemBC.hpp"

class EBC_Partition_2D : public IEBC_Partition
{
  public:
    EBC_Partition_2D( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const IElemBC * const &ebc_list );

    virtual ~EBC_Partition_2D();

    virtual void write_hdf5(const char * FileName) const;

    virtual void print_info() const;

  private:
    const int cpu_rank;

    std::vector<int> LFront_Elem, LBack_Elem, LLeft_Elem, LRight_Elem, LTop_Elem, LBottom_Elem;

    std::vector<int> Num_LBCElem;
};


#endif
