#ifndef IEN_FEM_HPP
#define IEN_FEM_HPP
// ============================================================================
// IEN_FEM.hpp
//
// This class defines the IEN array for a FEM mesh using uniform element type.
//
// According to the uniform element type assumption, nLocBas here is determined
// by in_ien.size() / in_nelem.
//
// Author: Ju Liu
// Date: July 30 2023
// ============================================================================
#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "IIEN.hpp"

class IEN_FEM final : public IIEN
{
  public:
    IEN_FEM( const int &in_nelem, const std::vector<int> &in_ien ) 
      : nElem( in_nelem ),
        nLocBas( VEC_T::get_size( in_ien )/nElem ),
        IEN( in_ien )
    {
      SYS_T::print_fatal_if( VEC_T::get_size( in_ien ) % nElem != 0,
          "Error: IEN_FEM input IEN vector size cannot do a perfect division by input number of elements. This suggests the mesh main contain different types of elements. \n" );
    }

    ~IEN_FEM() override = default;

    int get_IEN( const int &ee, const int &ii ) const override
    { return IEN[ee*nLocBas + ii]; }

    int get_nLocBas( const int &ee = 0 ) const override {return nLocBas;}

    void print_info() const override
    {
      std::cout<<std::endl;
      std::cout<<"====== IEN ====== \n";
      for(int ii=0; ii<nElem; ++ii)
      {
        for(int jj=0; jj<nLocBas; ++jj)
          std::cout<<get_IEN(ii, jj)<<'\t';
        std::cout<<std::endl;
      }
      std::cout<<"================= \n";
    }

  private:
    const int nElem, nLocBas;
    const std::vector<int> IEN;
};

#endif
