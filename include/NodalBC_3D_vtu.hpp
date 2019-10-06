#ifndef NODALBC_3D_VTU_HPP
#define NODALBC_3D_VTU_HPP
// ==================================================================
// NodalBC_3D_vtu.hpp
//
// This is an instantiation of INodalBC for 3D problems by reading a
// vtu file.
//
// This class is designed for fixing a group of nodes in a volumetric
// subdomain, e.g. the solid subdomain displacement in the FSI problem
// when solving the fluid mesh motion.
//
// Date: July 26 2017
// Author: Ju Liu
// ==================================================================
#include "INodalBC.hpp"
#include "Tet_Tools.hpp"

class NodalBC_3D_vtu : public INodalBC
{
  public:
    // --------------------------------------------------------------
    // Read in the vtu file that all the nodes in this file will be
    // enforced as essential boundary conditions.
    // --------------------------------------------------------------
    NodalBC_3D_vtu( const std::string &vtufilename, const int &nFunc );
   
    // --------------------------------------------------------------
    // Read in the vtu file and a list of vtp file that all nodes in these
    // files will be enforced as essential boundary conditions.
    // It is OK that the vtu and the vtp's have overlapping nodes.
    // --------------------------------------------------------------
    NodalBC_3D_vtu( const std::string &vtufilename, 
        const std::vector<std::string> &vtpfileList, const int &nFunc );
 
    // --------------------------------------------------------------
    // Read in the vtu file and two lists of vtp files.
    // The nodes in the vtu file extract the nodes
    // in the second vtp file, together with the nodes in the
    // vtpfileList will be the essential boundary conditions.
    // This is used for calculating the preloaded tissue state, while
    // the vtu file is the fluid domain, the vtpfileList are the inlet
    // and outlet faces, and the file to be extracted is the inner wall.
    // --------------------------------------------------------------
    NodalBC_3D_vtu( const std::string &vtufilename, 
        const std::vector<std::string> &vtpfileList,
        const std::vector<std::string> &vtpfileList_tobeextracted, 
        const int &nFunc );
 
    virtual ~NodalBC_3D_vtu();

  private:
    NodalBC_3D_vtu() {};
};

#endif
