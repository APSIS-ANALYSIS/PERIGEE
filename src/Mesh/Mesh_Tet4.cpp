#include "Mesh_Tet4.hpp"

Mesh_Tet4::Mesh_Tet4(const int &in_nfunc, const int &in_nelem)
  : nFunc(in_nfunc), nElem(in_nelem)
{}

Mesh_Tet4::~Mesh_Tet4()
{}

void Mesh_Tet4::print_info() const
{
  std::cout<<"======= Mesh_Tet4 ======="<<std::endl;
  std::cout<<"Degree: "<<get_s_degree()<<'\t'<<get_t_degree()<<'\t'<<get_u_degree()<<std::endl;
  std::cout<<"Total Elem: "<<get_nElem()<<std::endl;
  std::cout<<"Total Func: "<<get_nFunc()<<std::endl;
  std::cout<<"Local Basis #: "<<get_nLocBas()<<std::endl;
  std::cout<<"========================="<<std::endl;
}

// EOF
