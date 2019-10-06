#ifndef AGLOBAL_MESH_INFO_FEM_3D_HPP
#define AGLOBAL_MESH_INFO_FEM_3D_HPP
// ==================================================================
// AGlobal_Mesh_Info_FEM_3D.hpp
//
// This is the instantiation of global mesh info class which is used
// for 3D FEM mesh.
//
// Date Created: Jan 21 2017
// ==================================================================
#include "IAGlobal_Mesh_Info.hpp"
#include "HDF5_Reader.hpp"

class AGlobal_Mesh_Info_FEM_3D : public IAGlobal_Mesh_Info
{
  public:
    AGlobal_Mesh_Info_FEM_3D( const std::string &fileBaseName,
        const int &cpu_rank );

    // Construct a global mesh info based on a mesh by enriching it
    // in each cell with additional bubble nodes
    // Input:  num_enrich_node : the number of enriched nodes
    // Output: nFunc = original nFunc + nElem * num_enrich_node
    //         nLocBas = nLocBas + num_enrich_node
    //         elemType = elemType + 10
    AGlobal_Mesh_Info_FEM_3D( const std::string &fileBaseName,
        const int &cpu_rank, const int &num_enrich_node );
    
    virtual ~AGlobal_Mesh_Info_FEM_3D();

    virtual int get_xdegree() const {return xdegree;}
    virtual int get_ydegree() const {return ydegree;}
    virtual int get_zdegree() const {return zdegree;}

    virtual double get_max_hx() const
    {SYS_T::print_fatal("Error: AGlobal_Mesh_Info_FEM_3D::get_max_hx is not implemented. \n"); return 0.0;}

    virtual double get_max_hy() const
    {SYS_T::print_fatal("Error: AGlobal_Mesh_Info_FEM_3D::get_max_hy is not implemented. \n"); return 0.0;}

    virtual double get_max_hz() const
    {SYS_T::print_fatal("Error: AGlobal_Mesh_Info_FEM_3D::get_max_hz is not implemented. \n"); return 0.0;}
    
    virtual int get_nElem() const {return nElem;}

    virtual int get_nFunc() const {return nFunc;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual int get_probDim() const {return probDim;}

    virtual int get_elemType() const {return elemType;}

    virtual void print() const;

  private:
    int xdegree, ydegree, zdegree;
    int nElem, nFunc, nLocBas, probDim, elemType;
};

#endif
