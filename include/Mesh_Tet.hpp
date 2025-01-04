#ifndef MESH_TET_HPP
#define MESH_TET_HPP
// ============================================================================
// Mesh_Tet.hpp
//
// This is the instantiation of the IMesh class for tetrahedral element,
// i.e. unstructured mesh.
//
// Author: Ju Liu
// Date: July 30 2023
// ============================================================================
#include "IMesh.hpp"

class Mesh_Tet : public IMesh
{
  public:
    Mesh_Tet(const int &in_nFunc, const int &in_nElem, const int &in_deg) 
      : nFunc(in_nFunc), nElem(in_nElem)
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

    ~Mesh_Tet() override = default;

    void print_info() const override
    {
      std::cout<<"======= Mesh_Tet ======="<<std::endl;
      std::cout<<"Total Elem: "   <<get_nElem()<<std::endl;
      std::cout<<"Total Func: "   <<get_nFunc()<<std::endl;
      std::cout<<"Local Basis #: "<<get_nLocBas()<<std::endl;
      std::cout<<"========================="<<std::endl;
    }

    int get_degree() const override {return deg;}

    int get_nFunc() const override {return nFunc;}

    int get_nElem() const override {return nElem;}

    int get_nLocBas() const override {return nLocBas;}

  private:
    const int nFunc, nElem;

    // ------------------------------------------------------------------------ 
    // In this class, nLocBas is determined by the element type, and since we
    // are restricted to tet element, the degree determines the value of it.
    // ------------------------------------------------------------------------ 
    int nLocBas;
};

#endif
