#ifndef MESH_FEM_HPP
#define MESH_FEM_HPP
// ==================================================================
// Mesh_FEM.hpp
//
// This is the instantiation of the IMesh class for general element.
//
// Date: Nov. 17 2017
// ==================================================================
#include "IMesh.hpp"

class Mesh_FEM : public IMesh
{
  public:
    // Constructor
    // Assumes the polyminal degree is uniform
    Mesh_FEM(const int &in_nfunc, const int &in_nelem, const int &in_nlocbas, 
        const int &in_deg) : nFunc(in_nfunc), nElem(in_nelem), 
    nLocBas(in_nlocbas) {}

    ~Mesh_FEM() override = default;

    void print_info() const override
    {
      std::cout<<'\n';
      std::cout<<"======= Mesh_FEM ======="<<std::endl;
      std::cout<<"Total Elem: "<<get_nElem()<<std::endl;
      std::cout<<"Total Func: "<<get_nFunc()<<std::endl;
      std::cout<<"Local Basis #: "<<get_nLocBas()<<std::endl;
      std::cout<<"========================="<<std::endl;
    }

    int get_nFunc() const override {return nFunc;}

    int get_nElem() const override {return nElem;}

    int get_nLocBas() const override {return nLocBas;}

  private:
    const int nFunc, nElem, nLocBas;
};

#endif
