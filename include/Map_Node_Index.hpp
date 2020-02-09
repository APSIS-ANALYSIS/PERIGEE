#ifndef MAP_NODE_INDEX_HPP
#define MAP_NODE_INDEX_HPP
// ==================================================================
// Map_Node_Index.hpp
// Object:
// Stores the new indices that are reordered based on node partition.
// Stores two int vectors:
// old_2_new( ): map the old global index to new global index;
// new_2_old( ): map the new global index back to old global index.
//
// Date Created: Oct 3 2013
// ==================================================================
#include "IGlobal_Part.hpp"
#include "Vec_Tools.hpp"

class Map_Node_Index
{
  public:
    Map_Node_Index( const class IGlobal_Part * const &gpart,
        const int &cpu_size, const int nFunc );
    
    virtual ~Map_Node_Index();

    // Map the natural node numbering to the new numbering based on
    // mesh partition
    virtual int get_old2new(const int &ii) const {return old_2_new[ii];}
    
    // Map the new numbering back to the old, natural numbering for nodes
    virtual int get_new2old(const int &ii) const {return new_2_old[ii];}
    
    virtual void print_info() const;

    // ! write_hdf5(): write the old_2_new and new_2_old mappings as h5 file
    virtual void write_hdf5( const char * const &fileName ) const;

  private:
    std::vector<int> old_2_new, new_2_old;
};

#endif
