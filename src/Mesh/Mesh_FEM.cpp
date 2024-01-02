#include "Mesh_FEM.hpp"

Mesh_FEM::Mesh_FEM(const int &in_nfunc, const int &in_nelem,
    const int &in_nlocbas, const int &in_deg)
: nFunc(in_nfunc), nElem(in_nelem), nLocBas(in_nlocbas),
  sdeg(in_deg), tdeg(in_deg), udeg(in_deg)
{}

void Mesh_FEM::print_info() const
{
  std::cout<<'\n';
  std::cout<<"======= Mesh_FEM ======="<<std::endl;
  std::cout<<"Degree: "<<get_s_degree()<<'\t'<<get_t_degree()<<'\t'<<get_u_degree()<<std::endl;
  std::cout<<"Total Elem: "<<get_nElem()<<std::endl;
  std::cout<<"Total Func: "<<get_nFunc()<<std::endl;
  std::cout<<"Local Basis #: "<<get_nLocBas()<<std::endl;
  std::cout<<"========================="<<std::endl;
}

// EOF
