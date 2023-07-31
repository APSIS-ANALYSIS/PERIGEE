#ifndef IEN_GMSH_HPP
#define IEN_GMSH_HPP
// ==================================================================
// IEN_Gmsh.hpp
//
// This objects gives the IEN array generated from Gmsh.
// 
// The nodal indices in the reference element should be found from 
// the Gmsh manual.
//
// Date: Nov. 16 2017
// ==================================================================
#include <vector>
#include "IIEN.hpp"

class IEN_Gmsh : public IIEN
{
  public:
    IEN_Gmsh( const int &in_nelem, const int &in_nlocbas,
        const std::vector<int> &in_ien );

    ~IEN_Gmsh();

    virtual int get_IEN( const int &ee, const int &l_node ) const;

    virtual void print_info() const;

  private:
    const int nElem, nLocBas;

    int * IEN;
};

#endif
