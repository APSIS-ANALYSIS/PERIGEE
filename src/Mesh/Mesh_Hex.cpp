#include "Mesh_Hex.hpp"

Mesh_Hex::Mesh_Hex(const int &in_nfunc, const int &in_nelem, 
    const int &in_deg) : nFunc(in_nfunc), 
  nElem(in_nelem), deg(in_deg)
{
  switch(deg)
  {
    case 1:
      nLocBas = 8;
      break;
    case 2:
      nLocBas = 27;
      break;
    default:
      SYS_T::print_fatal("Error: Mesh_Hex, the input value of degree %d is not supported.\n", deg);
      nLocBas = 0;
      break;
  }
}

Mesh_Hex::~Mesh_Hex()
{}

void Mesh_Hex::print_info() const
{
  std::cout<<"======= Mesh_Hex ======="<<std::endl;
  std::cout<<"Degree: "       <<get_degree()<<std::endl;
  std::cout<<"Total Elem: "   <<get_nElem()<<std::endl;
  std::cout<<"Total Func: "   <<get_nFunc()<<std::endl;
  std::cout<<"Local Basis #: "<<get_nLocBas()<<std::endl;
  std::cout<<"========================="<<std::endl;
}

// EOF