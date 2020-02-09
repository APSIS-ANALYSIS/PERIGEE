#ifndef BC_PARTITION2D_HPP
#define BC_PARTITION2D_HPP
// ==================================================================
// BC_Partition2D.hpp
// Object:
// Partition the boundary nodes and boundary elements, and save them
// to .h5 files.
//
// Date:
// April 11 2014
// ==================================================================
#include "Vec_Tools.hpp"
#include "IPart.hpp"
#include "BoundaryCond2D.hpp"
#include "Map_Node_Index.hpp"

class BC_Partition2D
{
  public:
    BC_Partition2D( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        std::vector<BoundaryCond2D *> const &bc_list );
    
    virtual ~BC_Partition2D();

    virtual int get_LID(int pos) const {return LID[pos];}
    virtual int get_LDN(int pos) const {return LDN[pos];}
    virtual int get_LPSN(int pos) const {return LPSN[pos];}
    virtual int get_LPMN(int pos) const {return LPMN[pos];}
    virtual int get_Num_LD(int pos) const {return Num_LD[pos];}
    virtual int get_Num_LP(int pos) const {return Num_LP[pos];}

    virtual void write_hdf5(const char * FileName) const;
  private:
    int cpu_rank;

    // Local ID array
    std::vector<int> LID;
    
    // Local Dirichlet Nodes
    std::vector<int> LDN;

    // Local Periodic Slave / Master Nodes
    std::vector<int> LPSN, LPMN;

    // Number of Local Dirichlet Nodes
    std::vector<int> Num_LD;

    // Number of Local Periodic Slave Nodes
    std::vector<int> Num_LP;

    // List of boundary elements
    std::vector<int> LFront_Elem, LBack_Elem, LLeft_Elem, LRight_Elem, LTop_Elem, LBottom_Elem;
    
    // Number of Local Boundary Condition Element
    std::vector<int> Num_LBCElem;
};
#endif
