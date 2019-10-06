#include "BC_Partition.hpp"

BC_Partition::BC_Partition( const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    std::vector<BoundaryCond *> const &bc_list )
{
  cpu_rank = part->get_cpu_rank();
  // Initialize
  int dof = (int) bc_list.size();
  Num_LD.clear(); Num_LP.clear(); Num_LBCElem.clear();
  LID.clear(); LDN.clear(); LPSN.clear(); LPMN.clear();
  LFront_Elem.clear(); LBack_Elem.clear(); LLeft_Elem.clear();
  LRight_Elem.clear(); LTop_Elem.clear(); LBottom_Elem.clear();
  
  Num_LD.resize(dof); Num_LP.resize(dof);
  Num_LBCElem.resize(dof*6);// 6 faces for each dof
  
  // Loop of bc_list[i] and store dirichlet nodes, periodic
  // slave nodes that belong to this partition.
  int node_index, node_num;
  for(int i=0; i<dof; ++i)
  {
    node_num = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_dir_nodes(); ++j)
    {
      node_index = bc_list[i]->get_dir_nodes(j);
      node_index = mnindex->get_old2new(node_index);

      if(part->isNodeInPart(node_index) == true)
      {
        LDN.push_back(node_index);
        node_num += 1;
      }
    }
    Num_LD[i] = node_num;
  }

  int node_index_m;
  for(int i=0; i<dof; ++i)
  {
    node_num = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_per_nodes(); ++j)
    {
      node_index   = bc_list[i]->get_per_slave_nodes(j);
      node_index_m = bc_list[i]->get_per_master_nodes(j);
      node_index   = mnindex->get_old2new(node_index);
      node_index_m = mnindex->get_old2new(node_index_m);
      
      if(part->isNodeInPart(node_index) == true)
      {
        LPSN.push_back(node_index);
        LPMN.push_back(node_index_m);
        node_num += 1;
      }
    }
    Num_LP[i] = node_num;
  }
 
  int totnode = part->get_nlocghonode();

  LID.resize(totnode * dof);
 
  // ID array for the local nodes, including ghost nodes 
  for(int i=0; i<dof; ++i)
  {
    for(int j=0; j<totnode; ++j)
    {
      s_int new_index = part->get_local_to_global(j);
      s_int old_index = mnindex->get_new2old(new_index);
      LID[i*totnode+j] = bc_list[i]->get_ID( old_index );
    }
  }

  // For non dirichlet bc nodes, update the LID number to new global
  // node index
  for(int i=0; i<dof*totnode; ++i)
  {
    if(LID[i] != -1)
    {
      LID[i] = mnindex->get_old2new(LID[i]);
    }
  }

  VEC_T::shrink2fit(LID); VEC_T::shrink2fit(LDN);
  VEC_T::shrink2fit(LPSN); VEC_T::shrink2fit(LPMN);

  // Start Partition of boundary elements
  int lfro_elem, lbac_elem, llef_elem, lrig_elem, ltop_elem, lbot_elem;
  int temp_elem, elemlocindex;
  for(int i=0; i<dof; ++i)
  {
    lfro_elem = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_front_elem(); ++j)
    {
      temp_elem = bc_list[i]->get_front_elem(j);
      elemlocindex = part->get_elemLocIndex(temp_elem);
      if(elemlocindex != -1)
      {
        LFront_Elem.push_back(elemlocindex);
        lfro_elem += 1;
      }
    }
    Num_LBCElem[6*i + 0] = lfro_elem;
    
    lbac_elem = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_back_elem(); ++j)
    {
      temp_elem = bc_list[i]->get_back_elem(j);
      elemlocindex = part->get_elemLocIndex(temp_elem);
      if(elemlocindex != -1)
      {
        LBack_Elem.push_back(elemlocindex);
        lbac_elem += 1;
      }
    }
    Num_LBCElem[6*i + 1] = lbac_elem;

    llef_elem = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_left_elem(); ++j)
    {
      temp_elem = bc_list[i]->get_left_elem(j);
      elemlocindex = part->get_elemLocIndex(temp_elem);
      if(elemlocindex != -1)
      {
        LLeft_Elem.push_back(elemlocindex);
        llef_elem += 1;
      }
    }
    Num_LBCElem[6*i + 2] = llef_elem;

    lrig_elem = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_right_elem(); ++j)
    {
      temp_elem = bc_list[i]->get_right_elem(j);
      elemlocindex = part->get_elemLocIndex(temp_elem);
      if(elemlocindex != -1)
      {
        LRight_Elem.push_back(elemlocindex);
        lrig_elem += 1;
      }
    }
    Num_LBCElem[6*i + 3] = lrig_elem;

    ltop_elem = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_top_elem(); ++j)
    {
      temp_elem = bc_list[i]->get_top_elem(j);
      elemlocindex = part->get_elemLocIndex(temp_elem);
      if(elemlocindex != -1)
      {
        LTop_Elem.push_back(elemlocindex);
        ltop_elem += 1;
      }
    }
    Num_LBCElem[6*i + 4] = ltop_elem;

    lbot_elem = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_bottom_elem(); ++j)
    {
      temp_elem = bc_list[i]->get_bottom_elem(j);
      elemlocindex = part->get_elemLocIndex(temp_elem);
      if(elemlocindex != -1)
      {
        LBottom_Elem.push_back(elemlocindex);
        lbot_elem += 1;
      }
    }
    Num_LBCElem[6*i + 5] = lbot_elem;
  } 
}


BC_Partition::BC_Partition( const IMeshPart * const &part,
    const Map_Node_Index * const &mnindex,
    std::vector<BoundaryCond *> const &bc_list )
{
  cpu_rank = part->get_cpu_rank();
  // Initialize
  int dof = (int) bc_list.size();
  Num_LD.clear(); Num_LP.clear(); Num_LBCElem.clear();
  LID.clear(); LDN.clear(); LPSN.clear(); LPMN.clear();
  LFront_Elem.clear(); LBack_Elem.clear(); LLeft_Elem.clear();
  LRight_Elem.clear(); LTop_Elem.clear(); LBottom_Elem.clear();
  
  Num_LD.resize(dof); Num_LP.resize(dof);
  Num_LBCElem.resize(dof*6);// 6 faces for each dof
 
  // Loop of bc_list[i] and store dirichlet nodes, periodic
  // slave nodes that belong to this partition.
  int node_index, node_num;
  for(int i=0; i<dof; ++i)
  {
    node_num = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_dir_nodes(); ++j)
    {
      node_index = bc_list[i]->get_dir_nodes(j);
      node_index = mnindex->get_old2new(node_index);

      if(part->isNodeInPart(node_index) == true)
      {
        LDN.push_back(node_index);
        node_num += 1;
      }
    }
    Num_LD[i] = node_num;
  }

  int node_index_m;
  for(int i=0; i<dof; ++i)
  {
    node_num = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_per_nodes(); ++j)
    {
      node_index   = bc_list[i]->get_per_slave_nodes(j);
      node_index_m = bc_list[i]->get_per_master_nodes(j);
      node_index   = mnindex->get_old2new(node_index);
      node_index_m = mnindex->get_old2new(node_index_m);
      
      if(part->isNodeInPart(node_index) == true)
      {
        LPSN.push_back(node_index);
        LPMN.push_back(node_index_m);
        node_num += 1;
      }
    }
    Num_LP[i] = node_num;
  }
 
  int totnode = part->get_nlocghonode();

  LID.resize(totnode * dof);
 
  // ID array for the local nodes, including ghost nodes 
  for(int i=0; i<dof; ++i)
  {
    for(int j=0; j<totnode; ++j)
    {
      s_int new_index = part->get_local_to_global(j);
      s_int old_index = mnindex->get_new2old(new_index);
      LID[i*totnode+j] = bc_list[i]->get_ID( old_index );
    }
  }

  // For non dirichlet bc nodes, update the LID number to new global
  // node index
  for(int i=0; i<dof*totnode; ++i)
  {
    if(LID[i] != -1)
    {
      LID[i] = mnindex->get_old2new(LID[i]);
    }
  }

  VEC_T::shrink2fit(LID); VEC_T::shrink2fit(LDN);
  VEC_T::shrink2fit(LPSN); VEC_T::shrink2fit(LPMN);

  // Start Partition of boundary elements
  int lfro_elem, lbac_elem, llef_elem, lrig_elem, ltop_elem, lbot_elem;
  int temp_elem, elemlocindex;
  for(int i=0; i<dof; ++i)
  {
    lfro_elem = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_front_elem(); ++j)
    {
      temp_elem = bc_list[i]->get_front_elem(j);
      elemlocindex = part->get_elemLocIndex(temp_elem);
      if(elemlocindex != -1)
      {
        LFront_Elem.push_back(elemlocindex);
        lfro_elem += 1;
      }
    }
    Num_LBCElem[6*i + 0] = lfro_elem;
    
    lbac_elem = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_back_elem(); ++j)
    {
      temp_elem = bc_list[i]->get_back_elem(j);
      elemlocindex = part->get_elemLocIndex(temp_elem);
      if(elemlocindex != -1)
      {
        LBack_Elem.push_back(elemlocindex);
        lbac_elem += 1;
      }
    }
    Num_LBCElem[6*i + 1] = lbac_elem;

    llef_elem = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_left_elem(); ++j)
    {
      temp_elem = bc_list[i]->get_left_elem(j);
      elemlocindex = part->get_elemLocIndex(temp_elem);
      if(elemlocindex != -1)
      {
        LLeft_Elem.push_back(elemlocindex);
        llef_elem += 1;
      }
    }
    Num_LBCElem[6*i + 2] = llef_elem;

    lrig_elem = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_right_elem(); ++j)
    {
      temp_elem = bc_list[i]->get_right_elem(j);
      elemlocindex = part->get_elemLocIndex(temp_elem);
      if(elemlocindex != -1)
      {
        LRight_Elem.push_back(elemlocindex);
        lrig_elem += 1;
      }
    }
    Num_LBCElem[6*i + 3] = lrig_elem;

    ltop_elem = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_top_elem(); ++j)
    {
      temp_elem = bc_list[i]->get_top_elem(j);
      elemlocindex = part->get_elemLocIndex(temp_elem);
      if(elemlocindex != -1)
      {
        LTop_Elem.push_back(elemlocindex);
        ltop_elem += 1;
      }
    }
    Num_LBCElem[6*i + 4] = ltop_elem;

    lbot_elem = 0;
    for(unsigned int j=0; j<bc_list[i]->get_num_bottom_elem(); ++j)
    {
      temp_elem = bc_list[i]->get_bottom_elem(j);
      elemlocindex = part->get_elemLocIndex(temp_elem);
      if(elemlocindex != -1)
      {
        LBottom_Elem.push_back(elemlocindex);
        lbot_elem += 1;
      }
    }
    Num_LBCElem[6*i + 5] = lbot_elem;
  } 
}


BC_Partition::~BC_Partition()
{}

void BC_Partition::write_hdf5(const char * FileName) const
{
  // handle name for this partition
  std::string fName(FileName);
  fName.append("_p");

  if( cpu_rank / 10 == 0 )
    fName.append("0000");
  else if( cpu_rank / 100 == 0 )
    fName.append("000");
  else if( cpu_rank / 1000 == 0 )
    fName.append("00");
  else if( cpu_rank / 10000 == 0 )
    fName.append("0");

  std::stringstream sstrm;
  sstrm<<cpu_rank;
  fName.append(sstrm.str());

  fName.append(".h5");
 
  // hdf5 data tag 
  hid_t file_id, group_id;

  // open existing file that have part info
  file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  // create group
  group_id = H5Gcreate(file_id, "/bc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // LID
  hid_t dataspace_id_LID, setid_LID;
  hsize_t dim_LID[1];
  dim_LID[0] = (s_int) LID.size();
  dataspace_id_LID = H5Screate_simple(1, dim_LID, NULL);
  setid_LID = H5Dcreate( group_id, "LID", H5T_NATIVE_INT, dataspace_id_LID,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_LID, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &LID[0] );
  H5Dclose( setid_LID );
  H5Sclose( dataspace_id_LID );
  
  // LDN
  if(LDN.size() > 0)
  {
    hid_t dataspace_id_LDN, setid_LDN;
    hsize_t dim_LDN[1];
    dim_LDN[0] = (s_int) LDN.size();
    dataspace_id_LDN = H5Screate_simple(1, dim_LDN, NULL);
    setid_LDN = H5Dcreate( group_id, "LDN", H5T_NATIVE_INT, dataspace_id_LDN,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( setid_LDN, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &LDN[0] );
    H5Dclose( setid_LDN );
    H5Sclose( dataspace_id_LDN );
  }

  // LPSN, LPMN
  if(LPSN.size() > 0)
  {
    hid_t dataspace_id_LPN, setid_LPSN, setid_LPMN;
    hsize_t dim_LPN[1];
    dim_LPN[0] = (s_int) LPSN.size();
    dataspace_id_LPN = H5Screate_simple(1, dim_LPN, NULL);
    setid_LPSN = H5Dcreate( group_id, "LPSN", H5T_NATIVE_INT, dataspace_id_LPN,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( setid_LPSN, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &LPSN[0] );
    setid_LPMN = H5Dcreate( group_id, "LPMN", H5T_NATIVE_INT, dataspace_id_LPN,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( setid_LPMN, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &LPMN[0] );
    H5Dclose( setid_LPMN );
    H5Dclose( setid_LPSN );
    H5Sclose( dataspace_id_LPN );
  }

  // Num_LD
  hid_t dataspace_id_Num_LD, setid_Num_LD;
  hsize_t dim_Num_LD[1];
  dim_Num_LD[0] = (int) Num_LD.size();
  dataspace_id_Num_LD = H5Screate_simple(1, dim_Num_LD, NULL);
  setid_Num_LD = H5Dcreate( group_id, "Num_LD", H5T_NATIVE_INT, dataspace_id_Num_LD,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_Num_LD, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Num_LD[0] );
  H5Dclose( setid_Num_LD );
  H5Sclose( dataspace_id_Num_LD ); 

  // Num_LP
  hid_t dataspace_id_Num_LP, setid_Num_LP;
  hsize_t dim_Num_LP[1];
  dim_Num_LP[0] = (int) Num_LP.size();
  dataspace_id_Num_LP = H5Screate_simple(1, dim_Num_LP, NULL);
  setid_Num_LP = H5Dcreate( group_id, "Num_LP", H5T_NATIVE_INT, dataspace_id_Num_LP,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_Num_LP, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Num_LP[0] );
  H5Dclose( setid_Num_LP );
  H5Sclose( dataspace_id_Num_LP );

  // LFront_Elem
  if( LFront_Elem.size() > 0 )
  {
    hid_t dataspace_id_LFront_Elem, setid_LFront_Elem;
    hsize_t dim_LFront_Elem[1];
    dim_LFront_Elem[0] = (int) LFront_Elem.size();
    dataspace_id_LFront_Elem = H5Screate_simple(1, dim_LFront_Elem, NULL);
    setid_LFront_Elem = H5Dcreate( group_id, "LFront_Elem", H5T_NATIVE_INT, 
        dataspace_id_LFront_Elem, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( setid_LFront_Elem, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
        H5P_DEFAULT, &LFront_Elem[0] );
    H5Dclose( setid_LFront_Elem );
    H5Sclose( dataspace_id_LFront_Elem );
  }

  // LBack_Elem
  if( LBack_Elem.size() > 0 )
  {
    hid_t dataspace_id_LBack_Elem, setid_LBack_Elem;
    hsize_t dim_LBack_Elem[1];
    dim_LBack_Elem[0] = (int) LBack_Elem.size();
    dataspace_id_LBack_Elem = H5Screate_simple(1, dim_LBack_Elem, NULL);
    setid_LBack_Elem = H5Dcreate( group_id, "LBack_Elem", H5T_NATIVE_INT, 
        dataspace_id_LBack_Elem, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( setid_LBack_Elem, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
        H5P_DEFAULT, &LBack_Elem[0] );
    H5Dclose( setid_LBack_Elem );
    H5Sclose( dataspace_id_LBack_Elem );
  }

  // LLeft_Elem
  if( LLeft_Elem.size() > 0 )
  {
    hid_t dataspace_id_LLeft_Elem, setid_LLeft_Elem;
    hsize_t dim_LLeft_Elem[1];
    dim_LLeft_Elem[0] = (int) LLeft_Elem.size();
    dataspace_id_LLeft_Elem = H5Screate_simple(1, dim_LLeft_Elem, NULL);
    setid_LLeft_Elem = H5Dcreate( group_id, "LLeft_Elem", H5T_NATIVE_INT, 
        dataspace_id_LLeft_Elem, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( setid_LLeft_Elem, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
        H5P_DEFAULT, &LLeft_Elem[0] );
    H5Dclose( setid_LLeft_Elem );
    H5Sclose( dataspace_id_LLeft_Elem );
  }

  // LRight_Elem
  if( LRight_Elem.size() > 0 )
  {
    hid_t dataspace_id_LRight_Elem, setid_LRight_Elem;
    hsize_t dim_LRight_Elem[1];
    dim_LRight_Elem[0] = (int) LRight_Elem.size();
    dataspace_id_LRight_Elem = H5Screate_simple(1, dim_LRight_Elem, NULL);
    setid_LRight_Elem = H5Dcreate( group_id, "LRight_Elem", H5T_NATIVE_INT, 
        dataspace_id_LRight_Elem, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( setid_LRight_Elem, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
        H5P_DEFAULT, &LRight_Elem[0] );
    H5Dclose( setid_LRight_Elem );
    H5Sclose( dataspace_id_LRight_Elem );
  }

  // LTop_Elem
  if( LTop_Elem.size() > 0 )
  {
    hid_t dataspace_id_LTop_Elem, setid_LTop_Elem;
    hsize_t dim_LTop_Elem[1];
    dim_LTop_Elem[0] = (int) LTop_Elem.size();
    dataspace_id_LTop_Elem = H5Screate_simple(1, dim_LTop_Elem, NULL);
    setid_LTop_Elem = H5Dcreate( group_id, "LTop_Elem", H5T_NATIVE_INT, 
        dataspace_id_LTop_Elem, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( setid_LTop_Elem, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
        H5P_DEFAULT, &LTop_Elem[0] );
    H5Dclose( setid_LTop_Elem );
    H5Sclose( dataspace_id_LTop_Elem );
  }

  // LBottom_Elem
  if( LBottom_Elem.size() > 0 )
  {
    hid_t dataspace_id_LBottom_Elem, setid_LBottom_Elem;
    hsize_t dim_LBottom_Elem[1];
    dim_LBottom_Elem[0] = (int) LBottom_Elem.size();
    dataspace_id_LBottom_Elem = H5Screate_simple(1, dim_LBottom_Elem, NULL);
    setid_LBottom_Elem = H5Dcreate( group_id, "LBottom_Elem", H5T_NATIVE_INT, 
        dataspace_id_LBottom_Elem, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( setid_LBottom_Elem, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
        H5P_DEFAULT, &LBottom_Elem[0] );
    H5Dclose( setid_LBottom_Elem );
    H5Sclose( dataspace_id_LBottom_Elem );
  }

  hid_t dataspace_id_Num_LBCElem, setid_Num_LBCElem;
  hsize_t dim_Num_LBCElem[1];
  dim_Num_LBCElem[0] = (int) Num_LBCElem.size();
  dataspace_id_Num_LBCElem = H5Screate_simple(1, dim_Num_LBCElem, NULL);
  setid_Num_LBCElem = H5Dcreate( group_id, "Num_LBCElem", H5T_NATIVE_INT, 
      dataspace_id_Num_LBCElem, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_Num_LBCElem, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &Num_LBCElem[0] );
  H5Dclose( setid_Num_LBCElem );
  H5Sclose( dataspace_id_Num_LBCElem );


  // close the group and file
  H5Gclose(group_id);
  H5Fclose(file_id);
}

// EOF
