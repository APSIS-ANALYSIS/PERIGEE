#ifndef ALOCAL_NBC_HPP
#define ALOCAL_NBC_HPP
// ============================================================================
// ALocal_NBC.hpp
//
// FEM-analysis-use Local subdomain's Nodal Boundary Conditions.
//
// Author: Ju Liu
// Date: Oct. 8th 2015
// ============================================================================
#include "HDF5_Reader.hpp"

class ALocal_NBC
{
  public:
    // ------------------------------------------------------------------------
    // Read the Nodal BC info from the part file under the group name given by 
    // gname (nbc by default).
    // This function is specifically designed to read the nbc for mesh motion 
    // equations in FSI problem, which are stored under the gname mesh_nbc.
    // ------------------------------------------------------------------------
    ALocal_NBC( const std::string &fileBaseName, 
        const int &cpu_rank, const std::string &gname="/nbc" );

    virtual ~ALocal_NBC() = default;

    // ------------------------------------------------------------------------
    // print information
    // ------------------------------------------------------------------------
    virtual void print_info() const;

    // ------------------------------------------------------------------------
    // get the Local ID value for the node's dof-th degree-of-freedom.
    // 0 <= dof_index < dof , 0 <= node < nlocghonode
    // ------------------------------------------------------------------------
    virtual int get_LID(const int &dof_index, const int &node) const
    {return LID[dof_index * nlocghonode + node];}

    // ------------------------------------------------------------------------
    // get the Local ID value for the node assuming there is a single dof
    // attached to this node.
    // 0 <= node < nlocghonode
    // ------------------------------------------------------------------------
    virtual int get_LID( const int &node ) const {return LID[node];}

    // ------------------------------------------------------------------------
    // get_dof_LID: the Local ID array's dof number
    //              LID.size = dofMat * nlocghonode
    // ------------------------------------------------------------------------
    virtual int get_dof_LID() const {return dof;} 
    
    // ------------------------------------------------------------------------
    // get global indices of the Dirichlet nodes in the local subdomain
    // ------------------------------------------------------------------------
    virtual int get_LDN(const int &dof_index, const int &node) const
    {return LDN[LD_offset[dof_index] + node];}

    // ------------------------------------------------------------------------
    // get global indices of the slave nodes in the local subdomain
    // ------------------------------------------------------------------------
    virtual int get_LPSN(const int &dof_index, const int &ii) const
    {return LPSN[LPS_offset[dof_index] + ii];}

    // ------------------------------------------------------------------------
    // get global indices of the master nodes corresponding to the
    // slave nodes in the local subdomain, i.e. from get_LPSN() 
    // ------------------------------------------------------------------------
    virtual int get_LPMN(const int &dof_index, const int &ii) const
    {return LPMN[LPS_offset[dof_index] + ii];}

    // ------------------------------------------------------------------------
    // get global indices of the master nodes in the local subdomain 
    // ------------------------------------------------------------------------
    virtual int get_LocalMaster(const int &dof_index, const int &ii) const
    {return LocalMaster[LPM_offset[dof_index] + ii];}

    // ------------------------------------------------------------------------
    // get global indices of the slave nodes corresponding to the
    // master nodes in the local subdomain, i.e. from get_LocalMaster()
    // ------------------------------------------------------------------------
    virtual int get_LocalMasterSlave(const int &dof_index, const int &ii) const
    {return LocalMasterSlave[LPM_offset[dof_index] + ii];}

    // ------------------------------------------------------------------------
    // get the number of Dirichlet nodes in the local subdomain
    // ------------------------------------------------------------------------
    virtual int get_Num_LD(const int &dof_index) const 
    {return Num_LD[dof_index];}
    
    // ------------------------------------------------------------------------
    // get the number of slave nodes in the local subdomain
    // ------------------------------------------------------------------------
    virtual int get_Num_LPS(const int &dof_index) const 
    {return Num_LPS[dof_index];}

    // ------------------------------------------------------------------------
    // get the number of master nodes in the local subdomain
    // ------------------------------------------------------------------------
    virtual int get_Num_LPM(const int &dof_index) const 
    {return Num_LPM[dof_index];}

    // ------------------------------------------------------------------------
    // Clean the memory for local master nodes and their slaves, and vice versa 
    // ------------------------------------------------------------------------
    virtual void clean_LocalMaster()
    {
      VEC_T::clean(LocalMaster); 
      VEC_T::clean(LocalMasterSlave);
      VEC_T::clean(Num_LPM); 
      VEC_T::clean(LPM_offset);
    }

  protected:
    // dof := LID.size() / nlocghonode
    int dof, nlocghonode;

    // LID stores the ID of all degrees-of-freedom. The ID value could be the
    // corresponding grid point's index (old manner) or the corresponding matrix
    // column index (new manner). For the actual definition of the LID array,
    // the users should refer to the corresponding NBC_Partition routine.
    std::vector<int> LID;
   
    // LDN stores the Dirichlet nodes' index owned by this CPU locally;
    // LPSN stores the slave nodes' indices owned by this CPU locally;
    // LPMN stores the corresponding master nodes (of LPSN) indices owned by
    // this CPU locally;
    // LocalMaster stores the master nodes' indices owned by this CPU locally;
    // LocalMasterSlave stores the corresponding slave nodes (of LocalMaster)
    // indices owned by this CPU locally.
    // The meaning of the index meantioned above could be the grid point id (old
    // manner) or the corresponding matrix column index (new manner). The
    // actual definition is in the NBC_Partition routine. 
    std::vector<int> LDN, LPSN, LPMN, LocalMaster, LocalMasterSlave;

    std::vector<int> Num_LD, Num_LPS, Num_LPM;

    std::vector<int> LD_offset, LPS_offset, LPM_offset;
    
    ALocal_NBC() = delete; 
};

#endif
