#ifndef IEN_TETRA_P2_HPP
#define IEN_TETRA_P2_HPP
// ==================================================================
// IEN_Tetra_P2.hpp
//
// This class defines the IEN array for a full 3D quadratic
// tetrahedral mesh.
//
// The number of local nodes is 10, and the lenght of the IEN array
// is 10 x nElem.
//
// Author: Ju Liu
// Date: Jan 4 2020
// ==================================================================
#include <vector>
#include "IIEN.hpp"

class IEN_Tetra_P2 : public IIEN
{
  public:
    IEN_Tetra_P2( const int &in_nelem, const std::vector<int> &in_ien );

    ~IEN_Tetra_P2();

    virtual int get_IEN( const int &ee, const int &l_node ) const;

    virtual void print_IEN() const;

  private:
    const int nElem;
    
    int * IEN;
};

#endif
