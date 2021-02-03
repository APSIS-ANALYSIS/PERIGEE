#ifndef IEN_TETRA_P1_HPP
#define IEN_TETRA_P1_HPP
// ==================================================================
// IEN_Tetra_P1.hpp
//
// This objects gives the IEN array for a full 3D linear tetrahedral
// mesh.
//
// This mesh is assumed to be linear tetrahedral elements. Hence, the 
// number of local nodes is 4. The length of the IEN array is 4 nElem.
//
// Author: Ju Liu
// Date: Dec.18 2016.
// ==================================================================
#include <vector>
#include "IIEN.hpp"

class IEN_Tetra_P1 : public IIEN
{
  public:
    // Constructor: assume the ien_array is read from the .vtu file
    IEN_Tetra_P1( const int &in_nelem, const std::vector<int> &in_ien );
    
    ~IEN_Tetra_P1();

    virtual int get_IEN( const int &ee, const int &l_node ) const;

    virtual void print_IEN() const;

  private:
    const int nElem;
    
    int * IEN;
};

#endif
