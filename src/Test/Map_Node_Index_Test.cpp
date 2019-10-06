#include "Map_Node_Index_Test.hpp"

void Map_Node_Index_Test(const class Map_Node_Index * const &mnindex,
    const s_int &nFunc )
{
  for(s_int ii=0; ii<nFunc; ++ii)
  {
    s_int a = mnindex->get_old2new(ii);
    s_int b = mnindex->get_new2old(a);
    assert(b == ii);
  }

  for(s_int ii=0; ii<nFunc; ++ii)
  {
    s_int a = mnindex->get_new2old(ii);
    s_int b = mnindex->get_old2new(a);
    assert(b == ii);
  }
  cout<<"Map_Node_Index_Test: PASSED!"<<endl;
}

// EOF
