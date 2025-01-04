#ifndef IMESH_HPP
#define IMESH_HPP
// ==================================================================
// IMesh.hpp
// The Mesh class storet the basic parameters of the global FEM mesh.
// 
// Date: Sept 23rd 2013.
// ==================================================================
#include "Sys_Tools.hpp"

class IMesh
{
  public:
    IMesh( const int &in_nFunc, const int &in_nElem, const int &in_nlocbas) 
      : nFunc(in_nFunc), nElem(in_nElem), nLocBas(in_nlocbas) {};
    
    virtual ~IMesh() = default;

    void print_info() const
    {
      std::cout<<"========= IMesh ========="<<std::endl;
      std::cout<<"Total Elem: "<<get_nElem()<<std::endl;
      std::cout<<"Total Func: "<<get_nFunc()<<std::endl;
      std::cout<<"Local Basis Functions: "<<get_nLocBas()<<std::endl;
      std::cout<<"========================="<<std::endl;
    }

    int get_nFunc() const {return nFunc;}

    int get_nElem() const {return nElem;}

    int get_nLocBas() const {return nLocBas;}

  private:
    const int nFunc, nElem, nLocBas;
};

#endif
