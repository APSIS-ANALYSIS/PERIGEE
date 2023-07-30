#ifndef IEN_FEM_HPP
#define IEN_FEM_HPP
// ============================================================================
// IEN_FEM.hpp
//
// This class defines the IEN array for a FEM mesh using uniform element type.
//
// According to the uniform element type assumption, nLocBas here is determined
// by in_ien.size() / in_nelem.
//
// Author: Ju Liu
// Date: July 30 2023
// ============================================================================
#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "IIEN.hpp"

class IEN_FEM : public IIEN
{
  public:
    IEN_FEM( const int &in_nelem, const std::vector<int> &in_ien );

    ~IEN_FEM();

    virtual int get_IEN( const int &ee, const int &l_node ) const;

    virtual int get_nLocBas( const int &ee = 0 ) const {return nLocBas;}

    virtual void print_info() const;

  private:
    const int nElem, nLocBas;
    
    int * IEN;
};

#endif
