#ifndef AGLOBAL_MESH_INFO_FEM_3D_HPP
#define AGLOBAL_MESH_INFO_FEM_3D_HPP
// ==================================================================
// AGlobal_Mesh_Info_FEM_3D.hpp
//
// This is the global mesh info class used for a 3D FEM mesh.
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

    // Construct global mesh info by enriching each cell with 
    // additional bubble nodes
    // Input:  num_enrich_node : the number of enriched nodes per cell
    // Output: nFunc = original nFunc + nElem * num_enrich_node
    //         nLocBas = nLocBas + num_enrich_node
    //         elemType = elemType + 10
    AGlobal_Mesh_Info_FEM_3D( const std::string &fileBaseName,
        const int &cpu_rank, const int &num_enrich_node );
    
    virtual ~AGlobal_Mesh_Info_FEM_3D() = default;

    // Get the polynomial degree for the discretization method. For
    // unstructured mesh, the three function should return the same
    // value since one cannot differentiate the xyz direction in, e.g.,
    // tet mesh.
    virtual int get_xdegree() const {return xdegree;}
    virtual int get_ydegree() const {return ydegree;}
    virtual int get_zdegree() const {return zdegree;}

    // Get the total number of element of the whole mesh. 
    virtual int get_nElem() const {return nElem;}

    // Get the total number of nodes of the whole mesh.
    virtual int get_nFunc() const {return nFunc;}

    // Get the number of local basis functions for the element.
    // Note: this implicitly implies that we use the same type of element
    //       for the mesh.
    virtual int get_nLocBas() const {return nLocBas;}

    // Get the dimension of the problem.
    virtual int get_probDim() const {return probDim;}

    // Get an integer that indicate the element type.
    virtual int get_elemType() const {return elemType;}

    virtual void print_info() const;

  private:
    int xdegree, ydegree, zdegree;
    int nElem, nFunc, nLocBas, probDim, elemType;
    
    AGlobal_Mesh_Info_FEM_3D() = delete;
};

#endif
