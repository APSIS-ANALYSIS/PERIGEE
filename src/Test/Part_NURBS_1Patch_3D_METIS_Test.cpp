#include "Part_NURBS_1Patch_3D_METIS_Test.hpp"

void Part_NURBS_1Patch_3D_METIS_LIEN_Test(
    const class IPart * const &part,
    const class Map_Node_Index * const &mnindex,
    const class IIEN * const &IEN )
{
  for( int e = 0; e<part->get_nlocalele(); ++e )
  {
    for( int i = 0; i<part->get_nLocBas(); ++i )
    {
      s_int e_global = part->get_elem_loc(e);
      s_int index_ien = IEN->get_IEN(e_global, i);
      s_int index_lien = part->get_LIEN(e, i);

      index_ien = mnindex->get_old2new(index_ien);

      index_lien = part->get_local_to_global(index_lien);

      if( index_ien != index_lien )
      {
        cerr<<"ERROR: LIEN and IEN are incompatible at "<<e<<'\t'<<i<<endl;
        exit(1);
      }
    }
  }
  cout<<"LIEN Test: PASSED!"<<endl;
}

void Part_NURBS_1Patch_3D_METIS_Node_Test(
    const class IPart * const &part,
    const class Map_Node_Index * const &mnindex)
{
  for(int node=0; node<part->get_nlocalnode(); ++node)
  {
    s_int new_node = part->get_node_loc(node);
    s_int old_node = part->get_node_loc_original(node);

    if(new_node != mnindex->get_old2new(old_node))
    {
      cerr<<"ERROR: node_loc and node_loc_original are incompatible! \n";
      exit(1);
    }
    
    if(old_node != mnindex->get_new2old(new_node))
    {
      cerr<<"ERROR: node_loc and node_loc_original are incompatible! \n";
      exit(1);
    }
  }
  cout<<"node_loc and node_loc_original Test: PASSED! \n";

  for(int node=0; node<part->get_nlocalnode(); ++node)
  {
    assert(part->get_local_to_global(node) == part->get_node_loc(node));
  }
  for(int node=0; node<part->get_nghostnode(); ++node)
    assert(part->get_local_to_global(node+part->get_nlocalnode()) == part->get_node_ghost(node));

  cout<<"local_to_global Test: PASSED! \n";

  for(int e = 0; e<part->get_nlocalele()-1; ++e)
    assert(part->get_elem_loc(e) < part->get_elem_loc(e+1));
  for(int n=0; n<part->get_nlocalnode()-1; ++n)
  {
    assert(part->get_node_loc(n) < part->get_node_loc(n+1));
    assert(part->get_node_loc_original(n) < part->get_node_loc_original(n+1));
  }
  for(int n=0; n<part->get_nghostnode()-1; ++n)
    assert(part->get_node_ghost(n) < part->get_node_ghost(n+1));
  cout<<"Node is in assending order: PASSED! \n";
}

void Part_NURBS_1Patch_3D_METIS_Mesh_Test(
    const class IPart * const &part, const class IMesh * const &mesh)
{
  assert( part->get_sDegree() == mesh->get_s_degree() );
  assert( part->get_tDegree() == mesh->get_t_degree() );
  assert( part->get_uDegree() == mesh->get_u_degree() );
  
  assert( part->get_nElem() == mesh->get_nElem() );
  assert( part->get_nElem_x() == mesh->get_nElem_x() );
  assert( part->get_nElem_y() == mesh->get_nElem_y() );
  assert( part->get_nElem_z() == mesh->get_nElem_z() );
  
  assert( part->get_nFunc() == mesh->get_nFunc() );
  assert( part->get_nFunc_x() == mesh->get_nFunc_x() );
  assert( part->get_nFunc_y() == mesh->get_nFunc_y() );
  assert( part->get_nFunc_z() == mesh->get_nFunc_z() );
  
  assert( part->get_nLocBas() == mesh->get_nLocBas() );
  
  assert( part->get_hx_max() == mesh->get_hx_max() );
  assert( part->get_hy_max() == mesh->get_hy_max() );
  assert( part->get_hz_max() == mesh->get_hz_max() );
  
  assert( part->get_hx_min() == mesh->get_hx_min() );
  assert( part->get_hy_min() == mesh->get_hy_min() );
  assert( part->get_hz_min() == mesh->get_hz_min() );


  for(int e=0; e<part->get_nlocalele(); ++e)
  {
    s_int e_global = part->get_elem_loc(e);
    assert( part->get_hx(e) == mesh->get_hx(e_global));
    assert( part->get_hy(e) == mesh->get_hy(e_global));
    assert( part->get_hz(e) == mesh->get_hz(e_global));
  }

  cout<<"Mesh Test: PASSED! \n";
}

void Part_NURBS_1Patch_3D_METIS_CtrlPts_Test(
    const class IPart * const &part,
    const class Map_Node_Index * const &mnindex,
    const std::vector<double> &ctrlPts )
{
  for(int node=0; node<part->get_nlocalnode(); ++node)
  {
    double cp_x = part->get_ctrlPts_x_loc(node);
    double cp_y = part->get_ctrlPts_y_loc(node);
    double cp_z = part->get_ctrlPts_z_loc(node);
    double cp_w = part->get_ctrlPts_w_loc(node);

    s_int node_global = part->get_local_to_global(node);
    node_global = mnindex->get_new2old(node_global);

    double cpp_x = ctrlPts[node_global * 4 + 0];
    double cpp_y = ctrlPts[node_global * 4 + 1];
    double cpp_z = ctrlPts[node_global * 4 + 2];
    double cpp_w = ctrlPts[node_global * 4 + 3];

    assert(cp_x == cpp_x);
    assert(cp_y == cpp_y);
    assert(cp_z == cpp_z);
    assert(cp_w == cpp_w);
  }
  cout<<"Control Points Test: PASSED! \n";
}

void Part_NURBS_1Patch_3D_METIS_Print_Node(
    const class IMesh * const &mesh,
    const class Map_Node_Index * const &mnindex,
    const class IGlobal_Part * const &gpart )
{
  s_int nFunc_x = mesh->get_nFunc_x();
  s_int nFunc_y = mesh->get_nFunc_y();
  s_int nFunc_z = mesh->get_nFunc_z();

  std::ofstream efile("node_index.txt", std::ofstream::out | std::ofstream::trunc);

  for( s_int zz=0; zz<nFunc_z; ++zz )
  {
    efile<<"Level: "<<zz<<endl;
    for( s_int yy=nFunc_y-1; yy>=0; --yy )
    {
      for( s_int xx=0; xx<nFunc_x; ++xx )
      {
        s_int index = zz * nFunc_x * nFunc_y + yy * nFunc_x + xx;
        if(index / 10 == 0)
          efile<<0<<0<<index;
        else if(index / 100 == 0)
          efile<<0<<index;
        else
          efile<<index;

        efile<<"(";

        int part_index = gpart->get_npart(index);
        index = mnindex->get_old2new(index);
        if(index / 10 == 0)
          efile<<0<<0<<index;
        else if(index / 100 == 0)
          efile<<0<<index;
        else
          efile<<index;
        efile<<")"<<part_index<<"  ";
      }
      efile<<'\n';
    }
    efile<<'\n';
  }

  efile.close();
}

void Part_NURBS_1Patch_3D_METIS_Print_Elem(
    const class IMesh * const &mesh,
    const class IGlobal_Part * const &gpart )
{
   s_int nElem_x = mesh->get_nElem_x();
   s_int nElem_y = mesh->get_nElem_y();
   s_int nElem_z = mesh->get_nElem_z();
   std::ofstream efile("elem_part.txt", std::ofstream::out | std::ofstream::trunc);
   for( s_int zz = 0; zz<nElem_z; ++zz)
   {
     efile<<"Level: "<<zz<<endl;
     for( s_int yy=nElem_y-1; yy>=0; --yy )
     {
       for( s_int xx=0; xx<nElem_x; ++xx)
       {
         s_int index = zz * nElem_x * nElem_y + yy * nElem_x + xx;
        if(index / 10 == 0)
          efile<<0<<0<<index;
        else if(index / 100 == 0)
          efile<<0<<index;
        else
          efile<<index;

        efile<<"(";
        int part_index = gpart->get_epart(index);
        if(part_index / 10 == 0)
          efile<<0<<0<<part_index;
        else if(part_index / 100 == 0)
          efile<<0<<part_index;
        else
          efile<<part_index;
        efile<<")  ";
       }
       efile<<'\n';
     }
     efile<<'\n';
   }
   efile.close();
}
