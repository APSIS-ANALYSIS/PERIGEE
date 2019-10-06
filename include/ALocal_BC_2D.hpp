#ifndef ALOCAL_BC_2D_HPP
#define ALOCAL_BC_2D_HPP
// ==================================================================
// ALocal_BC_2D.hpp
// Local boundary condition for 2D analysis code use.
//
// Date:
// April 15 2014
// ==================================================================
#include "IALocal_BC.hpp"
#include "HDF5_PartReader.hpp"

#include "petscsys.h"

class ALocal_BC_2D : public IALocal_BC
{
  public:
    ALocal_BC_2D( const HDF5_PartReader * const &h5reader );
  
    virtual ~ALocal_BC_2D();
    
    virtual void print() const;

    // Returns the ID number for the degree-of-freedom dof_index
    // at node.
    virtual int get_LID(const int &dof_index, const int &node)
      const {return LID[dof_index * nlocghonode + node];}

    // Returns the LDN array for the degree-of-freedom dof_index
    // at the ii-th local dirichlet node.
    // ii ranges from 0 to Num_LD[dof_index]
    virtual int get_LDN(const int &dof_index, const int &ii) const;
    
    // Returns the LPSN array for the degree-of-freedom dof_index
    // at the ii-th local periodic slave nodes.
    // ii ranges from 0 to Num_LP[dof_index]
    virtual int get_LPSN(const int &dof_index, const int &ii) const;

    // Returns the local periodic master node
    virtual int get_LPMN(const int &dof_index, const int &ii) const;

    // Returns the local number of dirichlet nodes for the
    // degree-of-freedom dof_index
    virtual int get_Num_LD(const int &dof_index) const {return Num_LD[dof_index];}

    // Returns the local number of periodic nodes for the
    // degree-of-freedom dof_index
    virtual int get_Num_LP(const int &dof_index) const {return Num_LP[dof_index];}

    // Returns the local number of boundary element for the 
    // degree-of-freedom dof_index,
    // with its front face on the boundary
    virtual int get_NumLE_Fro(const int &dof_index) 
      const {return Num_LBCElem[6*dof_index + 0];}
  
    // with its back face on the boundary
    virtual int get_NumLE_Bac(const int &dof_index)
      const {return Num_LBCElem[6*dof_index + 1];}

    // with its left face on the boundary
    virtual int get_NumLE_Lef(const int &dof_index)
      const {return Num_LBCElem[6*dof_index + 2];}
    
    // with its right face on the boundary
    virtual int get_NumLE_Rig(const int &dof_index)
      const {return Num_LBCElem[6*dof_index + 3];}

    // with its left face on the boundary
    virtual int get_NumLE_Top(const int &dof_index)
      const {return Num_LBCElem[6*dof_index + 4];}

    // with its bottom face on the boundary
    virtual int get_NumLE_Bot(const int &dof_index)
      const {return Num_LBCElem[6*dof_index + 5];}

    // Returns the ii's local element index for the dof_index.
    // ii ranges from 0 to Num_LBCElem_XXX(dof_index)
    virtual int get_LFront_Elem(const int &dof_index, const int &ii) const;
    virtual int get_LBack_Elem(const int &dof_index, const int &ii) const;
    virtual int get_LLeft_Elem(const int &dof_index, const int &ii) const;
    virtual int get_LRight_Elem(const int &dof_index, const int &ii) const;
    virtual int get_LTop_Elem(const int &dof_index, const int &ii) const;
    virtual int get_LBottom_Elem(const int &dof_index, const int &ii) const;

  private:
    int dof, nlocghonode;
    std::vector<int> LID, LDN, LPSN, LPMN;
    std::vector<int> LFront_Elem, LBack_Elem, LLeft_Elem,
      LRight_Elem, LTop_Elem, LBottom_Elem, Num_LBCElem,
      Num_LD, Num_LP;
};
#endif
