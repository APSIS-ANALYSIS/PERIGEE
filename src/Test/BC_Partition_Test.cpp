#include "BC_Partition_Test.hpp"

void BC_Partition_Node_LID_Test( const BC_Partition * const &bcpart,
   const IPart * const &part, const int dof )
{
  int totnode = part->get_nlocghonode();
  for(int ii=0; ii<dof; ++ii)
  {
    for(int jj=0; jj<part->get_nlocalnode(); ++jj)
    {
      s_int g_index = part->get_local_to_global(jj);
      s_int i_index = bcpart->get_LID(ii*totnode+jj);
      if( i_index == -1 )
      {
        bool flag_d = isIndexInLDN(g_index, ii, bcpart);
        if(flag_d != true)
        {
          std::cerr<<"ERROR: "<<g_index<<" "<<i_index<<" at "<<ii<<" "<<jj<<std::endl;
          exit(1);
        }
      
      }
      else if( g_index != i_index )
      {
        bool flag_p = isIndexInLPSN(g_index, ii, bcpart);
        
        if(flag_p != true)
        {
          std::cerr<<"ERROR: "<<g_index<<" "<<bcpart->get_LID(ii*totnode+jj)<<" at "<<ii<<" "<<jj<<std::endl;
          exit(1);
        }

        flag_p = isIndexInLPMN(i_index, ii, bcpart);
        assert(flag_p == true);
      }
    }
  }
  std::cout<<"BC_Partition_Node_LID_Test: PASSED! \n";
}

void BC_Partition_Node_LID_Test_2( const BC_Partition * const &bcpart,
        const std::vector<BoundaryCond *> &bc_list, const IPart * const &part,
        const Map_Node_Index * const &mnindex, const int dof )
{
  int totnode = part->get_nlocghonode();
  for(int ii=0; ii<dof; ++ii)
  {
    BoundaryCond * bc = bc_list[ii];
    for(int jj=0; jj<totnode; ++jj)
    {
      s_int g_index = part->get_local_to_global(jj);
      s_int i_index = bcpart->get_LID(ii*totnode+jj);
      if( i_index == -1 )
      {
        bool flag_d = isIndexInDN(g_index, ii, bc, mnindex);
        if(flag_d != true)
        {
          std::cerr<<"ERROR: at "<<ii<<" "<<jj<<std::endl;
          exit(1);
        }
      }
      else if( g_index != i_index )
      {
        int pos;
        bool flag_p = isIndexInPSN(g_index, ii, bc, mnindex, pos);
        
        if(flag_p != true)
        {
          std::cerr<<"ERROR: at "<<ii<<" "<<jj<<std::endl;
          exit(1);
        }

        flag_p = isIndexInPMN(i_index, ii, bc, mnindex);
        assert(flag_p == true);
      }
    }
  }
  std::cout<<"BC_Partition_Node_LID_Test_2: PASSED! \n";
}

bool isIndexInLDN( const s_int index, const int dof, const BC_Partition * const &bcpart )
{
  int start = 0;
  for(int i=0; i<dof; ++i)
    start += bcpart->get_Num_LD(i);
  
  int eendd = start + bcpart->get_Num_LD(dof);

  bool flag = false;
  for(int ii= start; ii<eendd; ++ii)
  {
    if(bcpart->get_LDN(ii) == index)
      flag = true;
  }
  return flag;
}


bool isIndexInDN( const s_int index, const int dof, const BoundaryCond * const &bc,
   const Map_Node_Index * const &mnindex )
{
  int d_size = (int) bc->get_num_dir_nodes();
  for(int ii=0; ii<d_size; ++ii)
  {
    int dir_node = bc->get_dir_nodes(ii);
    if(index == mnindex->get_old2new(dir_node))
      return true;
  }
  return false;
}

bool isIndexInLPSN( const s_int index, const int dof, const BC_Partition * const &bcpart )
{
  int start = 0;
  for(int i=0; i<dof; ++i)
    start += bcpart->get_Num_LP(i);

  int eendd = start + bcpart->get_Num_LP(dof);
  bool flag = false;
  for(int ii= start; ii<eendd; ++ii)
  {
    if(bcpart->get_LPSN(ii) == index)
      flag = true;
  }
  return flag;
}

bool isIndexInPSN( const s_int index, const int dof, const BoundaryCond * const &bc,
   const Map_Node_Index * const &mnindex, int &pos )
{
  int d_size = (int) bc->get_num_per_nodes();
  for(int ii=0; ii<d_size; ++ii)
  {
    int per_node = bc->get_per_slave_nodes(ii);
    if(index == mnindex->get_old2new(per_node))
    {
      pos = ii;
      return true;
    }
  }
  return false;
}
  
bool isIndexInLPMN( const s_int index, const int dof,  const BC_Partition * const &bcpart )
{
  int start = 0;
  for(int i=0; i<dof; ++i)
    start += bcpart->get_Num_LP(i);

  int eendd = start + bcpart->get_Num_LP(dof);
  bool flag = false;
  for(int ii= start; ii<eendd; ++ii)
  {
    if(bcpart->get_LPMN(ii) == index)
      flag = true;
  }
  return flag;
}

bool isIndexInPMN( const s_int index, const int dof, const BoundaryCond * const &bc,
   const Map_Node_Index * const &mnindex )
{
  int d_size = (int) bc->get_num_per_nodes();
  for(int ii=0; ii<d_size; ++ii)
  {
    int per_node = bc->get_per_master_nodes(ii);
    if(index == mnindex->get_old2new(per_node))
      return true;
  }
  return false;
}

void BC_Partition_LID_Node_Test( const BC_Partition * const &bcpart,
    const std::vector<BoundaryCond *> &bc_list, const IPart * const &part,
    const Map_Node_Index * const &mnindex, const int dof )
{
  int totnode = part->get_nlocghonode();
  for(int ii=0; ii<dof; ++ii)
  {
    BoundaryCond * bc = bc_list[ii];
    for(int jj=0; jj<totnode; ++jj)
    {
      int g_index = part->get_local_to_global(jj);
      bool flag_d = isIndexInDN(g_index, ii, bc, mnindex);
      if(flag_d == true)
      {
        int i_index = bcpart->get_LID(ii*totnode+jj);
        assert(i_index == -1);
      }
      int pos;
      bool flag_p = isIndexInPSN(g_index, ii, bc, mnindex, pos);
      if(flag_p == true)
      {
        int i_index = bcpart->get_LID(ii*totnode+jj);
        int master = bc->get_per_master_nodes(pos);
        master = mnindex->get_old2new(master);
        assert(master == i_index);
      }
    }
  }
  std::cout<<"BC_Partition_LID_Node_Test: Passed! \n";
}


// EOF
