#ifndef IALOCAL_BC_HPP
#define IALOCAL_BC_HPP
// ==================================================================
// IALocal_BC.hpp
// Interface for Analysis Local boundary condition.
//
// Date:
// Nov. 10 2013
// ==================================================================
#include <iostream>
class IALocal_BC
{
  public:
    IALocal_BC(){};
  
    virtual ~IALocal_BC(){};

    virtual void print() 
      const {std::cout<<"This print() is not implemented. \n";}

    virtual int get_LID(const int &dof_index, const int &node) const = 0;

    virtual int get_LDN(const int &dof_index, const int &ii) const = 0;
    
    virtual int get_LPSN(const int &dof_index, const int &ii) const = 0;

    virtual int get_LPMN(const int &dof_index, const int &ii) const = 0;

    virtual int get_Num_LD(const int &dof_index) const = 0;

    virtual int get_Num_LP(const int &dof_index) const = 0;

    // Boundary element for 3-dimensional problems
    virtual int get_NumLE_Fro(const int &dof_index) const {return -1;} 
    virtual int get_NumLE_Bac(const int &dof_index) const {return -1;}
    virtual int get_NumLE_Lef(const int &dof_index) const {return -1;}
    virtual int get_NumLE_Rig(const int &dof_index) const {return -1;}
    virtual int get_NumLE_Top(const int &dof_index) const {return -1;}
    virtual int get_NumLE_Bot(const int &dof_index) const {return -1;}

    virtual int get_LFront_Elem(const int &dof_index, const int &ii) const {return -1;}
    virtual int get_LBack_Elem(const int &dof_index, const int &ii) const {return -1;}
    virtual int get_LLeft_Elem(const int &dof_index, const int &ii) const {return -1;}
    virtual int get_LRight_Elem(const int &dof_index, const int &ii) const {return -1;}
    virtual int get_LTop_Elem(const int &dof_index, const int &ii) const {return -1;}
    virtual int get_LBottom_Elem(const int &dof_index, const int &ii) const {return -1;}


};
#endif
