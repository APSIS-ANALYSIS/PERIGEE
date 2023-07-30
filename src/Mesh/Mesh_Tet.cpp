#include "Mesh_Tet10.hpp"

Mesh_Tet10::Mesh_Tet10(const int &in_nfunc, const int &in_nelem)
  : nFunc(in_nfunc), nElem(in_nelem)
{}

Mesh_Tet10::~Mesh_Tet10()
{}

void Mesh_Tet10::print_info() const
{
  std::cout<<"======= Mesh_Tet10 ======="<<std::endl;
  std::cout<<"Degree: "<<get_s_degree()<<'\t'<<get_t_degree()<<'\t'<<get_u_degree()<<std::endl;
  std::cout<<"Total Elem: "<<get_nElem()<<std::endl;
  std::cout<<"Total Func: "<<get_nFunc()<<std::endl;
  std::cout<<"Local Basis #: "<<get_nLocBas()<<std::endl;
  std::cout<<"========================="<<std::endl;
}

// EOF
