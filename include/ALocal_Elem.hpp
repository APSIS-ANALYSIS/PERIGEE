#ifndef ALOCAL_ELEM_HPP
#define ALOCAL_ELEM_HPP
// ==================================================================
// ALocal_Elem.hpp
// Analysis-use Local Element class. This class records the local
// partition's total number of elements and their global indices.
//
// Author: Ju Liu
// Date: Nov. 10 2013
// ==================================================================
#include "HDF5_Reader.hpp"

class ALocal_Elem
{
  public:
    // Constructor : read the h5 file of the given base name and rank
    ALocal_Elem(const std::string &fbasename, const int &cpu_rank);

    // Destructor
    virtual ~ALocal_Elem();

    // Return the element index based on the local element index.
    // 0 <= index < nlocalele
    virtual int get_elem_loc(const int &index) const {return elem_loc[index];}
   
    // Return the number of elements in this sub-domain owned by this CPU. 
    virtual int get_nlocalele() const {return nlocalele;}

    // Given the global element index, return its location in the vector
    // elem_loc. If it does not belong to this sub-domain, it will return -1
    virtual int get_pos(const int &global_e_index) const 
    {return VEC_T::get_pos(elem_loc, global_e_index);}

    virtual void print_info() const;

    // This is a virtual function for multiphysics simulations. A tag
    // is attached to each element to denote different physical domains,
    // such as fluid vs. solid subdomains. For a single domain problem,
    // this function is NOT needed, and returns a default value of 0.
    virtual int get_elem_tag(const int &index) const
    {
      SYS_T::print_fatal("Warning: ALocal_Elem::get_elem_tag is not implemented.\n");
      return 0;
    }

  private:
    // Global indices of elements that belong to the local CPU.
    std::vector<int> elem_loc;
    
    // The number of elements that belong to the CPU, which equals
    // the length of the elem_loc vector.
    int nlocalele;
};

#endif
