#ifndef IALOCAL_IEN_HPP
#define IALOCAL_IEN_HPP

// ============================================================================
// IALocal_IEN.hpp
//
// Interface for analysis-used Local IEN array.
// 
// This class is a pure virtual interface that provides a way to obtain LIEN
// arrays.
//
// Date: Aug. 24 2015
// ============================================================================
#include "Sys_Tools.hpp"

class IALocal_IEN
{
  public:
    IALocal_IEN(){};

    virtual ~IALocal_IEN(){};

    // ------------------------------------------------------------------------
    // get_LIEN: get the element ee's IEN array and pass it as a vector
    // ------------------------------------------------------------------------
    virtual void get_LIEN(const int &ee, std::vector<int> &elem_ien) const
    {SYS_T::commPrint("Error: get_LIEN is not implemented. \n");}

    
    // ------------------------------------------------------------------------
    // get_LIEN: get the element ee's IEN array and apss it as a dynamic array.
    // Note: users are responsible for allocating and deleting the dynamic
    // array.
    // ------------------------------------------------------------------------
    virtual void get_LIEN(const int &ee, int * &elem_ien) const
    {SYS_T::commPrint("Error: get_LIEN is not implemented. \n");}

    
    // ------------------------------------------------------------------------
    // get_nLocBas: the element ee's number of basis functions.
    //              0 <= ee < nlocalele = ALocal_Elem->get_nlocalele.
    // ------------------------------------------------------------------------
    virtual int get_nLocBas( const int &ee ) const
    {SYS_T::commPrint("Error: get_nLocBas is not implemented. \n"); return 0;}


    // ------------------------------------------------------------------------
    // isNode_in_Elem returns a bool that determines if the node belongs to the
    // element elem. The node index has to be the node index for the local
    // subdomain.
    // ------------------------------------------------------------------------
    virtual bool isNode_in_Elem(const int &elem, const int &node) const
    {SYS_T::commPrint("Error: isNode_in_Elem is not implemented. \n"); return false;}


    virtual void print_info() const
    {SYS_T::commPrint("Error: print_info is not implemented. \n");}

};



#endif
