#include "ALocal_BC_3D.hpp"

ALocal_BC_3D::ALocal_BC_3D(const HDF5_PartReader * const &h5reader)
{
  h5reader->get_BC_LID_dof(LID, dof);
  h5reader->get_BC_LD(LDN, Num_LD);
  h5reader->get_BC_LP(LPSN, LPMN, Num_LP);
  h5reader->get_BC_LBCE(Num_LBCElem, LFront_Elem, LBack_Elem,
      LLeft_Elem, LRight_Elem, LTop_Elem, LBottom_Elem);

  int a, b, c, d;
  std::vector<int> e, f, g, h;

  h5reader->get_LN(a, b, c, nlocghonode, d, e, f, g, h);
}

ALocal_BC_3D::~ALocal_BC_3D()
{}

int ALocal_BC_3D::get_LDN(const int &dof_index, const int &node) const
{
  int offset = 0;
  for(int ii=0; ii<dof_index; ++ii)
    offset += Num_LD[ii];

  return LDN[offset + node];
}

int ALocal_BC_3D::get_LPSN(const int &dof_index, const int &node) const
{
  int offset = 0;
  for(int ii=0; ii<dof_index; ++ii)
    offset += Num_LP[ii];

  return LPSN[offset + node];
}

int ALocal_BC_3D::get_LPMN(const int &dof_index, const int &node) const
{
  int offset = 0;
  for(int ii=0; ii<dof_index; ++ii)
    offset += Num_LP[ii];

  return LPMN[offset + node];
}

int ALocal_BC_3D::get_LFront_Elem(const int &dof_index, const int &ii) const
{
  int offset = 0;
  for(int jj=0; jj<dof_index; ++jj)
    offset += Num_LBCElem[6*jj + 0];
  return LFront_Elem[offset + ii];
}

int ALocal_BC_3D::get_LBack_Elem(const int &dof_index, const int &ii) const
{
  int offset = 0;
  for(int jj=0; jj<dof_index; ++jj)
    offset += Num_LBCElem[6*jj + 1];
  return LBack_Elem[offset + ii];
}

int ALocal_BC_3D::get_LLeft_Elem(const int &dof_index, const int &ii) const
{
  int offset = 0;
  for(int jj=0; jj<dof_index; ++jj)
    offset += Num_LBCElem[6*jj + 2];
  return LLeft_Elem[offset + ii];
}

int ALocal_BC_3D::get_LRight_Elem(const int &dof_index, const int &ii) const 
{
  int offset = 0;
  for(int jj=0; jj<dof_index; ++jj)
    offset += Num_LBCElem[6*jj + 3];
  return LRight_Elem[offset + ii];
}

int ALocal_BC_3D::get_LTop_Elem(const int &dof_index, const int &ii) const
{
  int offset = 0;
  for(int jj=0; jj<dof_index; ++jj)
    offset += Num_LBCElem[6*jj+4];
  return LTop_Elem[offset + ii];
}

int ALocal_BC_3D::get_LBottom_Elem(const int &dof_index, const int &ii) const
{
  int offset = 0;
  for(int jj=0; jj<dof_index; ++jj)
    offset += Num_LBCElem[6*jj+5];
  return LBottom_Elem[offset + ii];
}

void ALocal_BC_3D::print() const
{
  std::cout<<"ALocal_BC_3D: \n";

  std::cout<<"LID: \n";

  for(int ii=0; ii<dof; ++ii)
  {
    for(int jj=0; jj<nlocghonode; ++jj)
      std::cout<<get_LID(ii,jj)<<'\t';
    std::cout<<"\n \n";
  }
 
  std::cout<<"LDN: \n";
  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof "<<ii<<" : \n";
    for(int jj=0; jj<Num_LD[ii]; ++jj)
      std::cout<<get_LDN(ii, jj)<<'\t';
    std::cout<<"\n \n";
  }
 
  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof: "<<ii<<'\t';
    std::cout<<"Num_LD: "<<get_Num_LD(ii)<<'\t';
    std::cout<<"Num_LP: "<<get_Num_LP(ii)<<'\n';
  }
  
  std::cout<<std::endl;

  std::cout<<"LPSN: \n";
  for(int ii=0; ii<dof; ++ii)
  {
    for(int jj=0; jj<Num_LP[ii]; ++jj)
      std::cout<<get_LPSN(ii, jj)<<'\t';
    std::cout<<"\n \n";
  }

  std::cout<<"LPMN: \n";
  for(int ii=0; ii<dof; ++ii)
  {
    for(int jj=0; jj<Num_LP[ii]; ++jj)
      std::cout<<get_LPMN(ii, jj)<<'\t';
    std::cout<<"\n \n";
  }

  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof: "<<ii<<'\t';
    std::cout<<"NumLE_Fro: "<<get_NumLE_Fro(ii)<<'\t';
    std::cout<<"NumLE_Bac: "<<get_NumLE_Bac(ii)<<'\t';
    std::cout<<"NumLE_Lef: "<<get_NumLE_Lef(ii)<<'\t';
    std::cout<<"NumLE_Rig: "<<get_NumLE_Rig(ii)<<'\t';
    std::cout<<"NumLE_Top: "<<get_NumLE_Top(ii)<<'\t';
    std::cout<<"NumLE_Bot: "<<get_NumLE_Bot(ii)<<'\n';
  }


  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof: "<<ii<<'\n';
    std::cout<<"Front:  ";
    for(int jj=0; jj<get_NumLE_Fro(ii); ++jj)
      std::cout<<get_LFront_Elem(ii, jj)<<'\t';
    std::cout<<'\n';
    std::cout<<"Back:  ";
    for(int jj=0; jj<get_NumLE_Bac(ii); ++jj)
      std::cout<<get_LBack_Elem(ii, jj)<<'\t';
    std::cout<<'\n';
    std::cout<<"Left:  ";
    for(int jj=0; jj<get_NumLE_Lef(ii); ++jj)
      std::cout<<get_LLeft_Elem(ii, jj)<<'\t';
    std::cout<<'\n';
    std::cout<<"Right:  ";
    for(int jj=0; jj<get_NumLE_Rig(ii); ++jj)
      std::cout<<get_LRight_Elem(ii, jj)<<'\t';
    std::cout<<'\n';
    std::cout<<"Top:  ";
    for(int jj=0; jj<get_NumLE_Top(ii); ++jj)
      std::cout<<get_LTop_Elem(ii, jj)<<'\t';
    std::cout<<'\n';
    std::cout<<"Bottom: ";
    for(int jj=0; jj<get_NumLE_Bot(ii); ++jj)
      std::cout<<get_LBottom_Elem(ii, jj)<<'\t';
    std::cout<<std::endl;
  }
}
