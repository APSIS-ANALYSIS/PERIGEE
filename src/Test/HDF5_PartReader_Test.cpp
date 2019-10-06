#include "HDF5_PartReader_Test.hpp"

void HDF5_PartReader_Print_GMI( const class HDF5_PartReader * const &preader )
{
  cout<<"Global Mesh Info: \n";

  int sDegree, tDegree, uDegree;
  preader->get_GMI_degree(sDegree, tDegree, uDegree);
  cout<<"degree: "<<sDegree<<'\t'<<tDegree<<'\t'<<uDegree<<endl;

  double hx_max, hy_max, hz_max;
  double hx_min, hy_min, hz_min;
  preader->get_GMI_h_max(hx_max, hy_max, hz_max);
  preader->get_GMI_h_min(hx_min, hy_min, hz_min);

  cout<<"h max: "<<hx_max<<'\t'<<hy_max<<'\t'<<hz_max<<endl;
  cout<<"h min: "<<hx_min<<'\t'<<hy_min<<'\t'<<hz_min<<endl;

  int nElem, nElem_x, nElem_y, nElem_z;
  int nFunc, nFunc_x, nFunc_y, nFunc_z;
  preader->get_GMI_nElem(nElem, nElem_x, nElem_y, nElem_z);
  preader->get_GMI_nFunc(nFunc, nFunc_x, nFunc_y, nFunc_z);

  cout<<"nElem: "<<nElem<<'\t'<<nElem_x<<'\t'<<nElem_y<<'\t'<<nElem_z<<endl;
  cout<<"nFunc: "<<nFunc<<'\t'<<nFunc_x<<'\t'<<nFunc_y<<'\t'<<nFunc_z<<endl;
  
  int nLocBas;
  preader->get_GMI_nLocBas( nLocBas );
  cout<<"nLocBas: "<<nLocBas<<endl;
}

void HDF5_PartReader_Print_EXT( const class HDF5_PartReader * const &preader,
   const int &e )
{
  cout<<"Local Element "<<e<<": "<<endl;
  std::vector<double> ext_x, ext_y, ext_z;

  //preader->get_EXT_x(e, ext_x);
  //preader->get_EXT_y(e, ext_y);
  //preader->get_EXT_z(e, ext_z);

  cout<<"x-dir: "<<endl;
  VEC_T::print(ext_x);

  cout<<"y-dir: "<<endl;
  VEC_T::print(ext_y);

  cout<<"z-dir: "<<endl;
  VEC_T::print(ext_z);
}

void HDF5_PartReader_Print_LIEN( const class HDF5_PartReader * const &preader,
    const int &e )
{
  cout<<"Local Element "<<e<<": "<<endl;
  std::vector<int> LIEN;
  preader->get_LIEN(e, LIEN);
  cout<<"len: "<<endl;
  VEC_T::print(LIEN);
}

void HDF5_PartReader_Print_LE( const class HDF5_PartReader * const &preader )
{
  cout<<"Local Element: \n";
  std::vector<int> elem_loc; int nlocalele;
  preader->get_LE(elem_loc, nlocalele);
  cout<<"elem_loc: "<<nlocalele<<endl;
  VEC_T::print(elem_loc);
}

void HDF5_PartReader_Print_LN( const class HDF5_PartReader * const &preader )
{
  int nlocalnode, nghostnode, nbadnode, ntotalnode, nlocghonode;
  std::vector<int> local_to_global, node_ghost, node_loc, node_loc_original;

  preader->get_LN( nlocalnode, nghostnode, nbadnode, nlocghonode, ntotalnode, local_to_global,
    node_ghost, node_loc, node_loc_original );

  cout<<"Local_Node: \n";
  cout<<"nlocalnode: "<<nlocalnode<<endl;
  cout<<"nghostnode: "<<nghostnode<<endl;
  cout<<"nbadnode: "<<nbadnode<<endl;
  cout<<"nlocghonode: "<<nlocghonode<<endl;
  cout<<"ntotalnode: "<<ntotalnode<<endl;

  cout<<"local_to_global: \n";
  VEC_T::print(local_to_global);
  cout<<"node_loc: \n";
  VEC_T::print(node_loc);
  cout<<"node_loc_original: \n";
  VEC_T::print(node_loc_original);
  cout<<"node_ghost: \n";
  VEC_T::print(node_ghost);
}

void HDF5_PartReader_Print_PI( const class HDF5_PartReader * const &preader )
{
  int cpu_size, cpu_rank, dual_edge_ncommon;
  preader->get_PI(cpu_rank, cpu_size, dual_edge_ncommon);
  cout<<"Part_Info: \n";
  cout<<"cpu_rank: "<<cpu_rank<<" cpu_size: "<<cpu_size<<" dual_edge_ncommon: "
    <<dual_edge_ncommon<<endl;
}

void HDF5_PartReader_Print_CPL( const class HDF5_PartReader * const &preader )
{
  std::vector<double> ctrl_x, ctrl_y, ctrl_z, ctrl_w;
  preader->get_CPL(ctrl_x, ctrl_y, ctrl_z, ctrl_w);
  cout<<"ctrlPts_x_loc: "<<endl;
  VEC_T::print(ctrl_x);
  cout<<"ctrlPts_y_loc: "<<endl;
  VEC_T::print(ctrl_y);
  cout<<"ctrlPts_z_loc: "<<endl;
  VEC_T::print(ctrl_z);
  cout<<"ctrlPts_w_loc: "<<endl;
  VEC_T::print(ctrl_w);
}


void HDF5_PartReader_Print_BC_LID_dof( const class HDF5_PartReader * const &preader )
{
  std::vector<int> LID;
  int dof;
  preader->get_BC_LID_dof(LID, dof);
  cout<<"dof: "<<dof<<endl;
  cout<<"LID: "<<endl;
  VEC_T::print(LID);
}

void HDF5_PartReader_Print_BC_LDN( const class HDF5_PartReader * const &preader )
{
  std::vector<int> LDN, Num_LD;
  preader->get_BC_LD(LDN, Num_LD);
  cout<<"Num_LD: ";
  VEC_T::print(Num_LD);
  cout<<'\n'<<"LDN: \n";
  VEC_T::print(LDN);
}

void HDF5_PartReader_Print_BC_LP( const class HDF5_PartReader * const &preader )
{
  std::vector<int> lpmn, lpsn, num_lp;
  preader->get_BC_LP(lpsn, lpmn, num_lp);
  cout<<"Num_LP: ";
  VEC_T::print(num_lp);
  cout<<"LPSN: "<<endl;
  VEC_T::print(lpsn);
  cout<<"LPMN: "<<endl;
  VEC_T::print(lpmn);
}

void HDF5_PartReader_Print_BC_BCE( const class HDF5_PartReader * const &preader )
{
  std::vector<int> num_lbce, front, back, left, right, top, bottom;
  preader->get_BC_LBCE(num_lbce, front, back, left, right, top, bottom);
  cout<<"Num_LBCElem: ";
  VEC_T::print(num_lbce);
  cout<<"Front: \n";
  VEC_T::print(front);
  cout<<"Back: \n";
  VEC_T::print(back);
  cout<<"Left: \n";
  VEC_T::print(left);
  cout<<"Right: \n";
  VEC_T::print(right);
  cout<<"Top: \n";
  VEC_T::print(top);
  cout<<"Bottom: \n";
  VEC_T::print(bottom);

}

void HDF5_PartReader_Check( const class HDF5_PartReader * const &preader, 
    const class IPart * const &part )
{
  HDF5_PartReader_Check_CtrlPts(preader, part);
  HDF5_PartReader_Check_LIEN(preader, part);
}

void HDF5_PartReader_Check_CtrlPts( const class HDF5_PartReader * const &preader, 
    const class IPart * const &part )
{
  std::vector<double> ctrl_x, ctrl_y, ctrl_z, ctrl_w;
  preader->get_CPL(ctrl_x, ctrl_y, ctrl_z, ctrl_w);
  
  for(unsigned int ii=0; ii<ctrl_x.size(); ++ii)
  {
    if(ctrl_x[ii] != part->get_ctrlPts_x_loc(ii))
    {
      cerr<<"ERROR: reader failed to read control point x. \n";
      exit(1);
    }
  }
  
  for(unsigned int ii=0; ii<ctrl_y.size(); ++ii)
  {
    if(ctrl_y[ii] != part->get_ctrlPts_y_loc(ii))
    {
      cerr<<"ERROR: reader failed to read control point y. \n";
      exit(1);
    }
  }

  for(unsigned int ii=0; ii<ctrl_z.size(); ++ii)
  {
    if(ctrl_z[ii] != part->get_ctrlPts_z_loc(ii))
    {
      cerr<<"ERROR: reader failed to read control point z. \n";
      exit(1);
    }
  }
  
  for(unsigned int ii=0; ii<ctrl_w.size(); ++ii)
  {
    if(ctrl_w[ii] != part->get_ctrlPts_w_loc(ii))
    {
      cerr<<"ERROR: reader failed to read control point w. \n";
      exit(1);
    }
  }

  cout<<"Control points have been read correctly. \n";
}

void HDF5_PartReader_Check_LIEN( const class HDF5_PartReader * const &preader, 
    const class IPart * const &part )
{
  int nlocalele;
  std::vector<int> elem_loc;
  preader->get_LE(elem_loc, nlocalele);
  
  int nLocBas;
  preader->get_GMI_nLocBas(nLocBas);

  for(int ee=0; ee<nlocalele; ++ee)
  {
    std::vector<int> LIEN;
    preader->get_LIEN(ee, LIEN);
    for(int ii=0; ii<nLocBas; ++ii)
    {
      if(LIEN[ii] != part->get_LIEN(ee, ii))
      {
        cerr<<"ERROR: failed to read LIEN. \n";
        exit(1);
      }
    }
  }
  cout<<"LIEN has been read correctly. \n";
}

void HDF5_PartReader_Check_LID( const class HDF5_PartReader * const &preader, 
    const class BC_Partition * const &bcpart )
{
  std::vector<int> LID;
  int dof;
  preader->get_BC_LID_dof(LID, dof);

  for( unsigned int ii=0; ii<LID.size(); ++ii)
  {
    if( LID[ii] != bcpart->get_LID(ii) )
    {
      cerr<<"ERROR: LID read failed. \n";
      exit(1);
    }
  }
  cout<<"LID has been read correctly. \n";
}















