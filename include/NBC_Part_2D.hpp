#ifndef NBC_PART_2D_HPP
#define NBC_PART_2D_HPP
// ============================================================================
// NBC_Part_2D.hpp
//
// Nodal Boundary condition partition for two-dimensional meshes.
//
// Date: Oct. 7th 2015
// ============================================================================

#include "INBC_Partition.hpp"
#include "IPart.hpp"
#include "Map_Node_Index.hpp"
#include "INodalBC.hpp"
#include "Vec_Tools.hpp"
#include "HDF5_Writer.hpp"

class NBC_Part_2D : public INBC_Partition
{
  public:
    NBC_Part_2D( IPart const * const &part,
        Map_Node_Index const * const &mnindex,
        const vector<INodalBC *> &nbc_list );

    virtual ~NBC_Part_2D();


    virtual void write_hdf5( char const * const &FileName ) const;


    virtual void print_info() const;

  private:
    // The ID array for the local nodes
    vector<int> LID;

    // Local subdomain's Dirichlet node indices
    vector<unsigned int> LDN;


    // Local subdomain's slave nodes and their master nodes
    vector<unsigned int> LPSN, LPMN;


    // Local subdomain's master nodes and their slave nodes
    vector<unsigned int> LocalMaster, LocalMasterSlave;


    // Number of local dirichlet nodes for each D.O.F.
    vector<unsigned int> Num_LD;


    // Number of periodic slave nodes in local subdomain for each D.O.F.
    vector<unsigned int> Num_LPS;


    // Number of periodic master nodes in local subdomain for each D.O.F.
    vector<unsigned int> Num_LPM;

};

#endif
