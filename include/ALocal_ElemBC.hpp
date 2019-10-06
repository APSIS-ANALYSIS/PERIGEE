#ifndef ALOCAL_ELEMBC_HPP
#define ALOCAL_ELEMBC_HPP
// ============================================================================
// ALocal_ElemBC.hpp
//
// Analysis-use Local subdomain's Elemental Boundary Conditions.
//
// Date: Oct. 12 2015
// ============================================================================
#include "HDF5_PartReader.hpp"

class ALocal_ElemBC
{
  public:
    ALocal_ElemBC( const HDF5_PartReader * const &h5r );

    virtual ~ALocal_ElemBC();

    virtual void print_info() const;

    // Get the number of element on the boundary faces
    virtual int get_NumLE_Fro() const {return Num_LBCElem[0];}

    virtual int get_NumLE_Bac() const {return Num_LBCElem[1];}

    virtual int get_NumLE_Lef() const {return Num_LBCElem[2];}

    virtual int get_NumLE_Rig() const {return Num_LBCElem[3];}

    virtual int get_NumLE_Top() const {return Num_LBCElem[4];}

    virtual int get_NumLE_Bot() const {return Num_LBCElem[5];}


    // Get the element indices
    virtual int get_LFront_Elem(const int &ii) const {return LFront_Elem[ii];}

    virtual int get_LBack_Elem(const int &ii) const {return LBack_Elem[ii];}

    virtual int get_LLeft_Elem(const int &ii) const {return LLeft_Elem[ii];}

    virtual int get_LRight_Elem(const int &ii) const {return LRight_Elem[ii];}

    virtual int get_LTop_Elem(const int &ii) const {return LTop_Elem[ii];}

    virtual int get_LBottom_Elem(const int &ii) const {return LBottom_Elem[ii];}


  private:
    std::vector<int> Num_LBCElem;

    std::vector<int> LFront_Elem, LBack_Elem, LLeft_Elem, LRight_Elem, 
      LTop_Elem, LBottom_Elem;
};

#endif
