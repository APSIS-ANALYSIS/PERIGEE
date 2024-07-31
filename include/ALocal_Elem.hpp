#ifndef ALOCAL_ELEM_HPP
#define ALOCAL_ELEM_HPP
// ============================================================================
// ALocal_Elem.hpp
// Analysis-use local element class. This class records the local partition's 
// total number of elements, their global indices, and additional info.
//
// Author: Ju Liu
// Date: Nov. 10 2013
// ============================================================================
#include "HDF5_Reader.hpp"

class ALocal_Elem
{
  public:
    // ------------------------------------------------------------------------
    // Constructor : read the h5 file of the given base name and rank
    // ------------------------------------------------------------------------
    ALocal_Elem(const std::string &fbasename, const int &cpu_rank);

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    virtual ~ALocal_Elem() = default;

    // ------------------------------------------------------------------------
    // Return the element index based on the local element index.
    // 0 <= index < nlocalele
    // ------------------------------------------------------------------------
    virtual int get_elem_loc(const int &index) const {return elem_loc[index];}
   
    // ------------------------------------------------------------------------
    // Return the number of elements in this sub-domain owned by this CPU. 
    // ------------------------------------------------------------------------
    virtual int get_nlocalele() const {return nlocalele;}

    // ------------------------------------------------------------------------
    // Return the number of elements with tag value being the input tag_val.
    // This function can only be called when isTagged = true.
    // ------------------------------------------------------------------------
    virtual int get_nlocalele( const int &tag_val ) const;

    // ------------------------------------------------------------------------
    // Given the global element index, return its location in the vector
    // elem_loc. If it does not belong to this sub-domain, it will return -1
    // ------------------------------------------------------------------------
    virtual int get_pos(const int &global_e_index) const 
    {return VEC_T::get_pos(elem_loc, global_e_index);}

    virtual void print_info() const;

    // ------------------------------------------------------------------------
    // This is a virtual function for multiphysics simulations. A tag
    // is attached to each element to denote different physical domains,
    // such as fluid vs. solid subdomains. For a single domain problem,
    // this function is NOT needed, and returns a default value of 0.
    // ------------------------------------------------------------------------
    virtual int get_elem_tag(const int &ee) const
    {
      // if( isTagged ) return elem_tag[ee];
      // else return 0;
	  ASSERT(isTagged, "Error: get_elem_tag function 'isTagged' is faluse.\n");
      return elem_tag[ee];
    }

  private:
    // ------------------------------------------------------------------------
    // The number of elements that belong to the CPU, which equals the length 
    // of the elem_loc vector.
    // ------------------------------------------------------------------------
    int nlocalele;

    // ------------------------------------------------------------------------
    // Global indices of elements that belong to the local CPU.
    // ------------------------------------------------------------------------
    std::vector<int> elem_loc {};
    
    // ------------------------------------------------------------------------
    // Flag that determine if the element has an additional tag
    // ------------------------------------------------------------------------
    bool isTagged;

    // ------------------------------------------------------------------------
    // A vector recording the tag of elements. Length is nlocalele.
    // In FSI problems, we assume tag 0 gives fluid element; 
    //                            tag 1 gives solid element.
    // elem_tag is cleared if isTagged = false
    // ------------------------------------------------------------------------
    std::vector<int> elem_tag {};

    // Disallow default constructor
    ALocal_Elem() = delete;
};

#endif
