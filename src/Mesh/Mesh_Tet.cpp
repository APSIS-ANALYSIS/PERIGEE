#include "Mesh_Tet.hpp"

Mesh_Tet::Mesh_Tet(const int &in_nfunc, const int &in_nelem, 
    const int &in_deg) : nFunc(in_nfunc), 
  nElem(in_nelem), deg(in_deg)
{
  switch(deg)
  {
    case 1:
      nLocBas = 4;
      break;
    case 2:
      nLocBas = 10;
      break;
    default:
      SYS_T::print_fatal("Error: Mesh_Tet, the input value of degree %d is not supported.\n", deg);
      nLocBas = 0;
      break;
  }
}

void Mesh_Tet::print_info() const
{
  std::cout<<"======= Mesh_Tet ======="<<std::endl;
  std::cout<<"Degree: "       <<get_degree()<<std::endl;
  std::cout<<"Total Elem: "   <<get_nElem()<<std::endl;
  std::cout<<"Total Func: "   <<get_nFunc()<<std::endl;
  std::cout<<"Local Basis #: "<<get_nLocBas()<<std::endl;
  std::cout<<"========================="<<std::endl;
}

// EOF
