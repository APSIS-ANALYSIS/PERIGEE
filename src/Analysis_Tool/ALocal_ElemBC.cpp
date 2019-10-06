#include "ALocal_ElemBC.hpp"

ALocal_ElemBC::ALocal_ElemBC( const HDF5_PartReader * const &h5r )
{
  h5r->get_EBC(Num_LBCElem, LFront_Elem, LBack_Elem, LLeft_Elem,
      LRight_Elem, LTop_Elem, LBottom_Elem);
}



ALocal_ElemBC::~ALocal_ElemBC()
{
  VEC_T::clean(Num_LBCElem);
  VEC_T::clean(LFront_Elem);
  VEC_T::clean(LBack_Elem);
  VEC_T::clean(LLeft_Elem);
  VEC_T::clean(LRight_Elem);
  VEC_T::clean(LTop_Elem);
  VEC_T::clean(LBottom_Elem);
}


void ALocal_ElemBC::print_info() const
{
  std::cout<<"NumLE_Fro: "<<get_NumLE_Fro()<<'\t';
  std::cout<<"NumLE_Bac: "<<get_NumLE_Bac()<<'\t';
  std::cout<<"NumLE_Lef: "<<get_NumLE_Lef()<<'\t';
  std::cout<<"NumLE_Rig: "<<get_NumLE_Rig()<<'\t';
  std::cout<<"NumLE_Top: "<<get_NumLE_Top()<<'\t';
  std::cout<<"NumLE_Bot: "<<get_NumLE_Bot()<<'\n';

  std::cout<<"Front:  ";
  for(int jj=0; jj<get_NumLE_Fro(); ++jj)
    std::cout<<get_LFront_Elem(jj)<<'\t';
  std::cout<<'\n';
  std::cout<<"Back:  ";
  for(int jj=0; jj<get_NumLE_Bac(); ++jj)
    std::cout<<get_LBack_Elem(jj)<<'\t';
  std::cout<<'\n';
  std::cout<<"Left:  ";
  for(int jj=0; jj<get_NumLE_Lef(); ++jj)
    std::cout<<get_LLeft_Elem(jj)<<'\t';
  std::cout<<'\n';
  std::cout<<"Right:  ";
  for(int jj=0; jj<get_NumLE_Rig(); ++jj)
    std::cout<<get_LRight_Elem(jj)<<'\t';
  std::cout<<'\n';
  std::cout<<"Top:  ";
  for(int jj=0; jj<get_NumLE_Top(); ++jj)
    std::cout<<get_LTop_Elem(jj)<<'\t';
  std::cout<<'\n';
  std::cout<<"Bottom: ";
  for(int jj=0; jj<get_NumLE_Bot(); ++jj)
    std::cout<<get_LBottom_Elem(jj)<<'\t';
  std::cout<<std::endl;
}

// EOF
