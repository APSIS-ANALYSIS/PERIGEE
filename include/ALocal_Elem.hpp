#ifndef ALOCAL_ELEM_HPP
#define ALOCAL_ELEM_HPP
// ==================================================================
// ALocal_Elem.hpp
// Analysis use Local Element class. This class records the local 
// partition domain's element index and total number of local elements.
//
// Author: Ju Liu
// Date: Nov. 10 2013
// ==================================================================
#include "HDF5_Reader.hpp"

class ALocal_Elem
{
  public:
    // Constructor : read h5 file by giving the part file base name and rank
    ALocal_Elem(const std::string &fileBaseName, const int &cpu_rank);

    virtual ~ALocal_Elem();

    virtual int get_elem_loc(const int &index) const {return elem_loc[index];}
    
    virtual int get_nlocalele() const {return nlocalele;}

    virtual void print_info() const;

    // This is a virtual function for multiphysics simulations. A tag
    // is attached to an element to identify different physical domain,
    // for example, fluid subdomain and solid subdomain.
    // For single domain problem, this function is NOT needed.
    virtual int get_elem_tag(const int &index) const
    {
      SYS_T::print_fatal("Error: ALocal_Elem::get_elem_tag is not implemented.\n");
      return -1;
    }

  private:
    // A vector storing the local cpu's element indices.
    std::vector<int> elem_loc;
    
    // The length of the elem_loc vector.
    int nlocalele;
};

#endif
