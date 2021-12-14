#ifndef MAP_NODE_INDEX_HPP
#define MAP_NODE_INDEX_HPP
// ==================================================================
// Map_Node_Index.hpp
// 
// Object:
// Stores the new indices that are reordered based on mesh partition.
// In particular, it stores two integer vectors:
// old_2_new( ) : map the old global index to new global index;
// new_2_old( ) : map the new global index to old global index.
// 
// old global index refers to the index read from the original mesh;
// new golbal index refers to the index realigned according to the 
// mesh partitioning.
//
// Author: Ju Liu
// Date Created: Oct 3 2013
// ==================================================================
#include "IGlobal_Part.hpp"
#include "HDF5_Writer.hpp"
#include "HDF5_Reader.hpp"

class Map_Node_Index
{
  public:
    // Construct the index mapping based on the partitioning
    Map_Node_Index( const IGlobal_Part * const &gpart,
        const int &cpu_size, const int &nFunc );

    // Read npart from node_start to node_end-1 and record the node
    // with cpu index equaling to cpu_rank.
    // e.g., n_start = 0, n_end = nFunc_p for pressure field
    //       n_start = nFunc_p, n_end = nFunc_p + nFunc_v for velo field
    Map_Node_Index( const IGlobal_Part * const &gpart,
        const int &cpu_size, const int &n_start, const int &n_end );

    // Load the index mapping from file on disk
    Map_Node_Index( const char * const &fileName );

    virtual ~Map_Node_Index();

    // Map the natural node numbering to the new numbering based on
    // mesh partition
    virtual int get_old2new(const int &ii) const {return old_2_new[ii];}
    
    // Map the new numbering back to the old, natural numbering for nodes
    virtual int get_new2old(const int &ii) const {return new_2_old[ii];}
    
    virtual void print_info() const;

    // write the old_2_new and new_2_old mappings into an hdf5 file
    virtual void write_hdf5( const char * const &fileName ) const;

  private:
    std::vector<int> old_2_new, new_2_old;
};

#endif
