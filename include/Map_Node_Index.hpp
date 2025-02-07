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
    // Construct the index mapping based on the partitioning for field 0
    Map_Node_Index( const IGlobal_Part * const &gpart,
        const int &cpu_size, const int &nFunc, const int &field = 0 );

    // Load the index mapping from file on disk
    Map_Node_Index( const char * const &fileName );

    virtual ~Map_Node_Index() = default;

    // Map the natural node numbering to the new numbering based on
    // mesh partition
    virtual int get_old2new(const int &ii) const {return old_2_new[ii];}
    
    // Map the new numbering back to the old, natural numbering for nodes
    virtual int get_new2old(const int &ii) const {return new_2_old[ii];}
    
    virtual void print_info() const;

    // write the old_2_new and new_2_old mappings into an hdf5 file
    virtual void write_hdf5( const std::string &fileName ) const;

  private:
    std::vector<int> old_2_new, new_2_old;
};

#endif
