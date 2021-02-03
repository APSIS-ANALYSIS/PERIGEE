#ifndef ALOCAL_NODALBC_HPP
#define ALOCAL_NODALBC_HPP
// ==================================================================
// ALocal_NodalBC.hpp
//
// FEM-analysis-use Local subdomain's Nodal Boundary Conditions.
//
// Author: Ju Liu
// Date: Oct. 8th 2015
// ===================================================================
#include "HDF5_Reader.hpp"

class ALocal_NodalBC
{
  public:
    // --------------------------------------------------------------
    // Read the Nodal BC info from the part file under the group
    // name given by gname.
    // This function is specifically designed to read the nbc for 
    // mesh motion equations in FSI problem, which are stored under
    // the name mesh_nbc.
    // --------------------------------------------------------------
    ALocal_NodalBC( const std::string &fileBaseName, 
        const int &cpu_rank, const std::string &gname="/nbc" );

    virtual ~ALocal_NodalBC();

    // ---------------------------------------------------------------
    // print information
    // ---------------------------------------------------------------
    virtual void print_info() const;

    virtual int get_LID(const int &dof_index, const int &node) const
    {return LID[dof_index * nlocghonode + node];}

    // ---------------------------------------------------------------
    // get dofMat: the implicit solver's dof number
    //             LID.size = dofMat * nlocghonode
    // ---------------------------------------------------------------
    virtual int get_dofMat() const {return dof;} 
    
    // ---------------------------------------------------------------
    // get global indices of the Dirichlet nodes in the local subdomain
    // ---------------------------------------------------------------
    virtual int get_LDN(const int &dof_index, const int &node) const
    {return LDN[LD_offset[dof_index] + node];}

    // ---------------------------------------------------------------
    // get global indices of the slave nodes in the local subdomain
    // ---------------------------------------------------------------
    virtual int get_LPSN(const int &dof_index, const int &ii) const
    {return LPSN[LPS_offset[dof_index] + ii];}

    // ---------------------------------------------------------------
    // get global indices of the master nodes corresponding to the
    // slave nodes in the local subdomain, i.e. from get_LPSN() 
    // ---------------------------------------------------------------
    virtual int get_LPMN(const int &dof_index, const int &ii) const
    {return LPMN[LPS_offset[dof_index] + ii];}

    // ---------------------------------------------------------------
    // get global indices of the master nodes in the local subdomain 
    // ---------------------------------------------------------------
    virtual int get_LocalMaster(const int &dof_index, const int &ii) const
    {return LocalMaster[LPM_offset[dof_index] + ii];}

    // ---------------------------------------------------------------
    // get global indices of the slave nodes corresponding to the
    // master nodes in the local subdomain, i.e. from get_LocalMaster()
    // ---------------------------------------------------------------
    virtual int get_LocalMasterSlave(const int &dof_index, const int &ii) const
    {return LocalMasterSlave[LPM_offset[dof_index] + ii];}

    // ---------------------------------------------------------------
    // get the number of Dirichlet nodes in the local subdomain
    // ---------------------------------------------------------------
    virtual int get_Num_LD(const int &dof_index) const 
    {return Num_LD[dof_index];}
    
    // ---------------------------------------------------------------
    // get the number of slave nodes in the local subdomain
    // ---------------------------------------------------------------
    virtual int get_Num_LPS(const int &dof_index) const 
    {return Num_LPS[dof_index];}

    // ---------------------------------------------------------------
    // get the number of master nodes in the local subdomain
    // ---------------------------------------------------------------
    virtual int get_Num_LPM(const int &dof_index) const 
    {return Num_LPM[dof_index];}

    // ---------------------------------------------------------------
    // Clean the memory for local master nodes and their slaves, and
    // vice versa 
    // ---------------------------------------------------------------
    virtual void clean_LocalMaster()
    {
      VEC_T::clean(LocalMaster); 
      VEC_T::clean(LocalMasterSlave);
      VEC_T::clean(Num_LPM); 
      VEC_T::clean(LPM_offset);
    }

  protected:
    int dof, nlocghonode;

    std::vector<int> LID, LDN, LPSN, LPMN, LocalMaster, LocalMasterSlave;

    std::vector<int> Num_LD, Num_LPS, Num_LPM;

    std::vector<int> LD_offset, LPS_offset, LPM_offset;
};

#endif
