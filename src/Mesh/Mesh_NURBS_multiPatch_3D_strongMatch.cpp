#include "Mesh_NURBS_multiPatch_3D_strongMatch.hpp"

Mesh_NURBS_multiPatch_3D_strongMatch::Mesh_NURBS_multiPatch_3D_strongMatch( 
    const std::vector<IMesh *> &in_mlist, const int &in_num_pat )
: numPat(in_num_pat), mlist(in_mlist)
{
  if(mlist.size() != numPat)
  {
    std::cerr<<"Error: number of patches != mlist size. \n";
    exit(EXIT_FAILURE);
  }
  
  hx_max = mlist[0]->get_hx_max();
  hx_min = mlist[0]->get_hx_min();

  hy_max = mlist[0]->get_hy_max();
  hy_min = mlist[0]->get_hy_min();

  hz_max = mlist[0]->get_hz_max();
  hz_min = mlist[0]->get_hz_min();

  for(unsigned int ii=0; ii<mlist.size(); ++ii)
  {
    if(hx_max<mlist[ii]->get_hx_max())
      hx_max = mlist[ii]->get_hx_max();
    if(hx_min>mlist[ii]->get_hx_min())
      hx_min = mlist[ii]->get_hx_min();
    
    if(hy_max<mlist[ii]->get_hy_max())
      hy_max = mlist[ii]->get_hy_max();
    if(hy_min>mlist[ii]->get_hy_min())
      hy_min = mlist[ii]->get_hy_min();
    
    if(hz_max<mlist[ii]->get_hz_max())
      hz_max = mlist[ii]->get_hz_max();
    if(hz_min>mlist[ii]->get_hz_min())
      hz_min = mlist[ii]->get_hz_min();
  }

  s_degree = mlist[0]->get_s_degree();
  t_degree = mlist[0]->get_t_degree();
  u_degree = mlist[0]->get_u_degree();
  nLocBas  = mlist[0]->get_nLocBas();

  for(unsigned int ii=0; ii<mlist.size(); ++ii)
  {
    if( s_degree != mlist[ii]->get_s_degree() )
    {
      std::cerr<<"Error: patch "<<ii<<" s_degree != "<<s_degree<<std::endl;
      exit(EXIT_FAILURE);
    }
    if( t_degree != mlist[ii]->get_t_degree() )
    {
      std::cerr<<"Error: patch "<<ii<<" t_degree != "<<t_degree<<std::endl;
      exit(EXIT_FAILURE);
    }
    if( u_degree != mlist[ii]->get_u_degree() )
    {
      std::cerr<<"Error: patch "<<ii<<" u_degree != "<<u_degree<<std::endl;
      exit(EXIT_FAILURE);
    }
    if( nLocBas != mlist[ii]->get_nLocBas() )
    {
      std::cerr<<"Error: patch "<<ii<<" nLocBas != "<<nLocBas<<std::endl;
      exit(EXIT_FAILURE);
    }
  }

  elem_ptr.clear();
  elem_ptr.resize(numPat+1);
  elem_ptr[0] = 0;

  nFunc = 0; nElem = 0;
  for(unsigned int ii = 0; ii<mlist.size(); ++ii)
  {
    nFunc += mlist[ii]->get_nFunc();
    nElem += mlist[ii]->get_nElem();

    elem_ptr[ii+1] = elem_ptr[ii] + mlist[ii]->get_nElem();
  }

  // Check the function and element start index for each patch
  int counter = 0;
  for(unsigned int ii=0; ii<mlist.size(); ++ii)
  {
    if(mlist[ii]->get_nFunc_start() != counter)
    {
      std::cerr<<"Error: Patch "<<ii<<" starting function index does not match";
      std::cerr<<" with the previous patches' nFunc. \n";
      exit(EXIT_FAILURE);
    }
    counter += mlist[ii]->get_nFunc();

    if(mlist[ii]->get_nElem_start() != elem_ptr[ii])
    {
      std::cerr<<"Error: Patch "<<ii<<" starting element index does not match";
      std::cerr<<" with the previous patches' nElem. \n";
      exit(EXIT_FAILURE);
    }
  }


  // Finish mesh format check
}



Mesh_NURBS_multiPatch_3D_strongMatch::~Mesh_NURBS_multiPatch_3D_strongMatch()
{
  std::cout<<"-- Mesh_NURBS_multiPatch_3D_strongMatch deleted. \n";
}



double Mesh_NURBS_multiPatch_3D_strongMatch::get_hx( s_int ee ) const
{
  int p_ind, e_ind;
  get_locelem_index(ee, p_ind, e_ind);
  return mlist[p_ind]->get_hx(e_ind);
}


double Mesh_NURBS_multiPatch_3D_strongMatch::get_hy( s_int ee ) const
{
  int p_ind, e_ind;
  get_locelem_index(ee, p_ind, e_ind);
  return mlist[p_ind]->get_hy(e_ind);
}



double Mesh_NURBS_multiPatch_3D_strongMatch::get_hz( s_int ee ) const
{
  int p_ind, e_ind;
  get_locelem_index(ee, p_ind, e_ind);
  return mlist[p_ind]->get_hz(e_ind);
}


void Mesh_NURBS_multiPatch_3D_strongMatch::get_locelem_index(const int &ee,
    int &pind, int &loc_ee) const
{
  pind = 0;
  while(ee>=elem_ptr[pind+1])
    ++pind;

  loc_ee = ee - elem_ptr[pind];
}



void Mesh_NURBS_multiPatch_3D_strongMatch::print_mesh_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== Mesh_NURBS_multiPatch_3D_strongMatch ======"<<std::endl;
  std::cout<<"Degree: "<<get_s_degree()<<'\t'<<get_t_degree()<<'\t'<<get_u_degree()<<std::endl;
  std::cout<<"Max mesh size: "<<get_hx_max()<<'\t'<<get_hy_max()<<'\t'<<get_hz_max()<<std::endl;
  std::cout<<"Min mesh size: "<<get_hx_min()<<'\t'<<get_hy_min()<<'\t'<<get_hz_min()<<std::endl;
  std::cout<<"Total Elem: "<<get_nElem()<<std::endl;
  std::cout<<"Total Func: "<<get_nFunc()<<std::endl;
  std::cout<<"Local Basis #: "<<get_nLocBas()<<std::endl;
  std::cout<<"elem_ptr : "<<std::endl;
  std::cout<<"number of patches : "<<get_num_patch()<<std::endl;
  VEC_T::print(elem_ptr);
  for(unsigned int ii=0; ii<numPat; ++ii)
  {
    std::cout<<"Patch "<<ii;
    get_patch_mesh(ii)->print_mesh_info();
  }
}

// EOF
