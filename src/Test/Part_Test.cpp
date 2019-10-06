#include "Part_Test.hpp"

void TEST_T::Part_LIEN_Test(
    const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN )
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


void TEST_T::Part_Node_Test(
    const IPart * const &part,
    const Map_Node_Index * const &mnindex)
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


void TEST_T::Part_CtrlPts_Test(
    const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const std::vector<double> &ctrlPts )
{
  for(int node=0; node<part->get_nlocalnode(); ++node)
  {
    double cp_x = part->get_ctrlPts_x_loc(node);
    double cp_y = part->get_ctrlPts_y_loc(node);
    double cp_z = part->get_ctrlPts_z_loc(node);

    s_int node_global = part->get_local_to_global(node);
    node_global = mnindex->get_new2old(node_global);

    double cpp_x = ctrlPts[node_global * 3 + 0];
    double cpp_y = ctrlPts[node_global * 3 + 1];
    double cpp_z = ctrlPts[node_global * 3 + 2];

    assert(cp_x == cpp_x);
    assert(cp_y == cpp_y);
    assert(cp_z == cpp_z);
  }
  cout<<"Control Points Test: PASSED! \n";
}


void TEST_T::Part_NBC_Test( const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const std::vector<INodalBC *> &nbc,
    const INBC_Partition * const &nbcpart,
    const int &dof )
{
  // Check the LID array
  const int totnode = part -> get_nlocghonode();
  
  // generate the offset for the arrays that store the nodal indices in a 1D
  // vector
  std::vector<int> LD_offset, LPSN_offset, LPMN_offset;

  LD_offset.resize(dof);   LD_offset[0] = 0;
  LPSN_offset.resize(dof); LPSN_offset[0] = 0;
  LPMN_offset.resize(dof); LPMN_offset[0] = 0;

  for(int ii=0; ii<dof-1; ++ii)
  {
    LD_offset[ii+1]   = LD_offset[ii] + nbcpart->get_Num_LD(ii);
    LPSN_offset[ii+1] = LPSN_offset[ii] + nbcpart->get_Num_LPS(ii);
    LPMN_offset[ii+1] = LPMN_offset[ii] + nbcpart->get_Num_LPM(ii);
  }

  for(int ii=0; ii<dof; ++ii)
  {
    std::vector<int> ldir, lsla, lmas;
    ldir.resize(nbcpart->get_Num_LD(ii));
    lsla.resize(nbcpart->get_Num_LPS(ii));
    lmas.resize(nbcpart->get_Num_LPS(ii));
    
    for(int jj=0; jj<nbcpart->get_Num_LD(ii); ++jj)
      ldir[jj] = nbcpart->get_LDN(jj + LD_offset[ii]);
    for(int jj=0; jj<nbcpart->get_Num_LPS(ii); ++jj)
    {
      lsla[jj] = nbcpart->get_LPSN(jj + LPSN_offset[ii]);
      lmas[jj] = nbcpart->get_LPMN(jj + LPSN_offset[ii]);
    }

    for(int jj=0; jj<part->get_nlocalnode(); ++jj)
    {
      int gindex = part->get_local_to_global(jj);
      int iindex = nbcpart->get_LID(ii*totnode+jj);

      if( iindex == -1 )
      {
        if( !VEC_T::is_invec(ldir, gindex) )
        {
          cerr<<"Error: Part_NBC_Test, Dir nodes at "<<ii<<" "<<jj<<endl;
          exit(EXIT_FAILURE);
        }
      }
      else if( gindex != iindex )
      {
        int iindex_pos = VEC_T::get_pos(lmas, iindex);
        int gindex_pos = VEC_T::get_pos(lsla, gindex);

        if( iindex_pos < 0 || gindex_pos < 0 )
        {
          cout<<iindex<<'\t'<<gindex<<'\t';
          cout<<iindex_pos<<'\t'<<gindex_pos<<endl;
          cerr<<"Error: Part_NBC_Test, slave and master nodes does not match.\n";
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  cout<<"LID test: PASSED! \n";

  for(int ii=0; ii<dof; ++ii)
  {
    std::vector<int> gdir, gsla, gmas;
    gdir.resize(nbc[ii]->get_num_dir_nodes());
    gsla.resize(nbc[ii]->get_num_per_nodes());
    gmas.resize(nbc[ii]->get_num_per_nodes());

    for(unsigned int jj=0; jj<nbc[ii]->get_num_dir_nodes(); ++jj) 
      gdir[jj] = mnindex->get_old2new( nbc[ii]->get_dir_nodes(jj) );
    for(unsigned int jj=0; jj<nbc[ii]->get_num_per_nodes(); ++jj)
    {
      gsla[jj] = mnindex->get_old2new( nbc[ii]->get_per_slave_nodes(jj) );
      gmas[jj] = mnindex->get_old2new( nbc[ii]->get_per_master_nodes(jj) );
    }
    
    std::vector<int> ldir, lsla, lmas, lmmm, lmms;
    ldir.resize(nbcpart->get_Num_LD(ii));
    lsla.resize(nbcpart->get_Num_LPS(ii));
    lmas.resize(nbcpart->get_Num_LPS(ii));
    lmmm.resize(nbcpart->get_Num_LPM(ii));
    lmms.resize(nbcpart->get_Num_LPM(ii));
    
    for(int jj=0; jj<nbcpart->get_Num_LD(ii); ++jj)
      ldir[jj] = nbcpart->get_LDN(jj + LD_offset[ii]);
    for(int jj=0; jj<nbcpart->get_Num_LPS(ii); ++jj)
    {
      lsla[jj] = nbcpart->get_LPSN(jj + LPSN_offset[ii]);
      lmas[jj] = nbcpart->get_LPMN(jj + LPSN_offset[ii]);
    }
    for(int jj=0; jj<nbcpart->get_Num_LPM(ii); ++jj)
    {
      lmmm[jj] = nbcpart->get_LocalMaster(jj + LPMN_offset[ii]);
      lmms[jj] = nbcpart->get_LocalMasterSlave(jj + LPMN_offset[ii]);
    }

    for(unsigned int jj=0; jj<ldir.size(); ++jj)
    {
      if( VEC_T::get_pos(gdir, ldir[jj]) < 0 )
      {
        cerr<<"Error: dof "<<ii<<" node "<<ldir[jj]
          <<" is not a Dirichlet node.\n";
        exit(EXIT_FAILURE);
      }
    }

    for( unsigned int jj=0; jj<lsla.size(); ++jj )
    {
      if( VEC_T::get_pos(gsla, lsla[jj]) < 0 )
      {
        cerr<<"Error: dof "<<ii<<" node "<<lsla[jj]
          <<"is not in global slave node list. \n";
        exit(EXIT_FAILURE);
      }

      if( lmas[jj] != gmas[VEC_T::get_pos(gsla, lsla[jj])] )
      {
        cerr<<"Error: local slave and its master are not found in global periodic node list. \n";
        exit(EXIT_FAILURE);
      }
    }

    for( unsigned int jj=0; jj<lmmm.size(); ++jj )
    {
      if( VEC_T::get_pos(gmas, lmmm[jj]) < 0 )
      {
        cerr<<"Error: dof "<<ii<<" node "<<lmmm[jj]
          <<"is not int global master node list. \n";
        exit(EXIT_FAILURE);
      }

      if( lmmm[jj] != gmas[VEC_T::get_pos(gsla, lmms[jj])] )
      {
        cout<<lmms[jj]<<'\t'<<gsla[VEC_T::get_pos(gmas, lmmm[jj])]<<'\n';
        cerr<<"Error: local master and its slave are not found in global periodic nodt list.\n";
        exit(EXIT_FAILURE);
      }
    }
  }

  cout<<"Boundary node test: PASSED! \n";
}


void TEST_T::EBCPart_AEBC_Test(
    const IEBC_Partition * const &ebc,
    const ALocal_EBC * const &ebc_part )
{
  assert(ebc->get_num_ebc() == ebc_part->get_num_ebc());
  
  const int num_ebc = ebc->get_num_ebc();

  for(int ii=0; ii<num_ebc; ++ii)
  {
    assert(ebc->get_num_local_node(ii) == ebc_part->get_num_local_node(ii));
    assert(ebc->get_num_local_cell(ii) == ebc_part->get_num_local_cell(ii));
    assert(ebc->get_cell_nLocBas(ii) == ebc_part->get_cell_nLocBas(ii));
    
    for(int jj=0; jj<ebc->get_num_local_node(ii); ++jj)
    {
      assert(ebc->get_local_pt_xyz(ii,jj*3) == ebc_part->get_local_pt_xyz(ii,jj*3));
      assert(ebc->get_local_pt_xyz(ii,jj*3+1) == ebc_part->get_local_pt_xyz(ii,jj*3+1));
      assert(ebc->get_local_pt_xyz(ii,jj*3+2) == ebc_part->get_local_pt_xyz(ii,jj*3+2));
      assert(ebc->get_local_global_node(ii,jj) == ebc_part->get_local_global_node(ii,jj));
      assert(ebc->get_local_node_pos(ii,jj) == ebc_part->get_local_node_pos(ii,jj));
    }

    for(int jj=0; jj<ebc->get_num_local_cell(ii); ++jj)
    {
      const int nlocbas = ebc->get_cell_nLocBas(ii);
      for(int kk=0; kk<nlocbas; ++kk)
      {
        assert(ebc->get_local_tri_ien(ii, nlocbas*jj+kk) == ebc_part->get_local_tri_ien(ii, nlocbas*jj+kk));
      }
      assert(ebc->get_local_global_cell(ii,jj) == ebc_part->get_local_global_cell(ii,jj));
    }
  }
  cout<<"ALocal_EBC verified! \n";
}


void TEST_T::Part_EBC_Test( const IEBC_Partition * const &ep,
    const IPart * const &part, const ElemBC * const &ebc )
{
  assert(ebc->get_num_ebc() == ep->get_num_ebc());

  const int num_ebc = ebc->get_num_ebc();

  for(int ii=0; ii<num_ebc; ++ii)
  {
    const int num_node = ep->get_num_local_node(ii);
    const int num_cell = ep->get_num_local_cell(ii);
    const int nlocbas = ep->get_cell_nLocBas(ii);
    
    assert(nlocbas == ebc->get_cell_nLocBas(ii));

    std::vector<int> gcell, gnode;
    gcell.clear(); gnode.clear();
    for(int jj=0; jj<ebc->get_num_cell(ii); ++jj)
      gcell.push_back(ebc->get_global_cell(ii,jj));
    for(int jj=0; jj<ebc->get_num_node(ii); ++jj)
      gnode.push_back(ebc->get_global_node(ii, jj));

    for(int cell=0; cell<num_cell; ++cell)
    {
      int cell_gindex = ep->get_local_global_cell(ii, cell);
      int cell_pos = VEC_T::get_pos(gcell, cell_gindex);
      // This cell's global index is found in the ebc list
      assert(cell_pos >= 0);
      // and is found in part
      assert( part->isElemInPart(cell_gindex) ); 
      
      // check IEN
      std::vector<int> lien_e;
      lien_e.resize(part->get_nLocBas());
      
      for(int ll=0; ll<part->get_nLocBas(); ++ll)
        lien_e[ll] = part->get_LIEN(part->get_elemLocIndex(cell_gindex), ll);

      std::vector<int> lien_tri; lien_tri.resize(nlocbas);
      for(int ll=0; ll<nlocbas; ++ll)
      {
        // index in the local nodal list
        int node = ep->get_local_tri_ien(ii, cell*nlocbas+ll);
        int enode = ebc->get_ien(ii, cell_pos, ll);
        
        // node's global volumetric index
        int node_gindex = ep->get_local_global_node(ii, node);
        int enode_gindex = ebc->get_global_node(ii, enode);

        lien_tri[ll] = ep->get_local_node_pos(ii, node);
        
        assert(node_gindex == enode_gindex);
      }
      
      int pos0 = VEC_T::get_pos(lien_e, lien_tri[0]);
      int pos1 = VEC_T::get_pos(lien_e, lien_tri[1]);
      int pos2 = VEC_T::get_pos(lien_e, lien_tri[2]);
      
      assert(pos0 >=0 && pos0 < part->get_nLocBas());
      assert(pos1 >=0 && pos1 < part->get_nLocBas());
      assert(pos2 >=0 && pos2 < part->get_nLocBas());
    }

    for(int nn=0; nn<num_node; ++nn)
    {
      // check xyz coordinates
      double nn_x = ep->get_local_pt_xyz(ii, nn*3);
      double nn_y = ep->get_local_pt_xyz(ii, nn*3+1);
      double nn_z = ep->get_local_pt_xyz(ii, nn*3+2);

      int node_pos = VEC_T::get_pos(gnode, ep->get_local_global_node(ii,nn));
      assert(node_pos >= 0);

      double gn_x = ebc->get_pt_xyz(ii, node_pos, 0);
      double gn_y = ebc->get_pt_xyz(ii, node_pos, 1);
      double gn_z = ebc->get_pt_xyz(ii, node_pos, 2);

      assert(nn_x == gn_x);
      assert(nn_y == gn_y);
      assert(nn_z == gn_z);
    }
  }

  cout<<"Elemental BC: PASSED! \n";
}


// EOF
