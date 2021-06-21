#ifndef ALOCAL_ELEM_WTAG_HPP
#define ALOCAL_ELEM_WTAG_HPP
// ============================================================================
// ALocal_Elem_wTag.hpp
//
// Analysis-use Local element class with physical tags to identify different 
// sub-domains, where the governing equations are different.
// 
// The additional data elem_tag is an integer array with length nlocalele. 
// elem_tag[ii] records the elem_loc[ii]'s physical domain tag.
//
// In FSI problems, we assume tag 0 gives fluid element; tag 1 gives solid element.
//
// Date: July 28 2017
// Author: Ju Liu
// ============================================================================
#include "ALocal_Elem.hpp"

class ALocal_Elem_wTag : public ALocal_Elem
{
  public:
    ALocal_Elem_wTag(const std::string &fileBaseName, const int &cpu_rank);

    virtual ~ALocal_Elem_wTag();

    // ------------------------------------------------------------------------
    // ! get_elem_tag : return the tag for the local element
    //   \para 0 <= ee < nlocalele
    // ------------------------------------------------------------------------
    virtual int get_elem_tag(const int &ee) const {return elem_tag[ee];}

    virtual void print_info() const;

  private:
    std::vector<int> elem_tag;
};

#endif
