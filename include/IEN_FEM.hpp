#ifndef IEN_FEM_HPP
#define IEN_FEM_HPP
// ==================================================================
// IEN_FEM.hpp
//
// This class defines the IEN array for a FEM mesh using uniform 
// element type.
//
// Author: Ju Liu
// Date: July 30 2023
// ==================================================================
#include <vector>
#include "IIEN.hpp"

class IEN_FEM : public IIEN
{
  public:
    IEN_FEM( const int &in_nelem, const std::vector<int> &in_ien );

    ~IEN_FEM();

    virtual int get_IEN( const int &ee, const int &l_node ) const;

    virtual void print_IEN() const;

  private:
    const int nElem;
    
    int * IEN;
};

#endif
