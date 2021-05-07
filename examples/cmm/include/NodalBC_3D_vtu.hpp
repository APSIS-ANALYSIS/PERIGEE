#ifndef NODALBC_3D_VTU_HPP
#define NODALBC_3D_VTU_HPP
// ============================================================================
// NodalBC_3D_vtu.hpp
//
// This is an instantiation of INodalBC for 3D problems by reading vtu
// files.
//
// This class is designed for (1) fixing a group of nodes in a volumetric
// subdomain, e.g. the solid subdomain displacement in the FSI problem
// when solving the fluid mesh motion (2) enforce essential bc that
// are from .vtu files, such as the quadratic triangle mesh.
//
// Date: July 26 2017
// Author: Ju Liu
// ============================================================================
#include "INodalBC.hpp"
#include "Tet_Tools.hpp"

class NodalBC_3D_vtu : public INodalBC
{
  public:
    // ------------------------------------------------------------------------
    // Default constructor: clear the dir_nodes, per_slave_nodes,
    // per_master_nodes; set num_dir_nodes = num_per_nodes = 0;
    // Set ID array based on the above `no-nodal' bc setting.
    // ------------------------------------------------------------------------
    NodalBC_3D_vtu( const int &nFunc );

    // ------------------------------------------------------------------------ 
    // Specify the Dirichlet nodes for CMM-type FSI simulations.
    // \para nbc_inflow : inflow nodes
    // \para nbc_ring   : ring nodes
    // \para comp       : the dof components ranges from 0 to 2 representing
    //                    x-, y-, and z-components.
    // ------------------------------------------------------------------------ 
    NodalBC_3D_vtu( const INodalBC * const &nbc_inflow,
        const INodalBC * const &nbc_ring,
        const int &comp, const int &nFunc );

    // ------------------------------------------------------------------------
    // Read in the vtu file that all the nodes in this file will be
    // enforced as essential boundary conditions.
    // ------------------------------------------------------------------------
    NodalBC_3D_vtu( const std::string &vtufilename, const int &nFunc );

    // ------------------------------------------------------------------------
    // Read a list of vtu files to specify the Dirichlet nodes.
    // No periodic BC nodes.
    // ------------------------------------------------------------------------
    NodalBC_3D_vtu( const std::vector<std::string> &vtufileList,
        const int &nFunc );

    virtual ~NodalBC_3D_vtu();

  private:
    NodalBC_3D_vtu() {};
};

#endif
