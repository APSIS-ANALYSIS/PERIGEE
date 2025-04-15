#ifndef NODALBC_HPP
#define NODALBC_HPP
// ==================================================================
// NodalBC.hpp
//
// This is an instantiation of INodalBC by reading the Dirichlet nodes 
// from a vtk file.
//
// The data contained in this class include:
// dir_nodes : the nodal indices for the Dirichlet nodes;
// num_dir_nodes : the number of the Dirichlet nodes, i.e., the length
//                 of the dir_nodes array;
// per_slave_nodes : the nodal indices for the slave nodes;
// per_master_nodes : the nodal indices for the master nodes;
// num_per_nodes : the number of periodic-type nodes, i.e., the length 
//                 of the per_slave_nodes / per_master_nodes.
//
// ID : the vector for the ID array, which is generated based on the 
//      nFunc, the number of total basis functions, and the dir_nodes
//      and per_slave_nodes/per_master_nodes. Once the dir_nodes, per_
//      xxx_nodes, and nFunc are given, the ID array will be generated
//      by the function create_ID.
//
// Date Created: Aug. 24 2023
// Author: Ju Liu
// ==================================================================
#include "INodalBC.hpp"
#include "Vec_Tools.hpp"
#include "VTK_Tools.hpp"

class NodalBC : public INodalBC
{
  public:
    // --------------------------------------------------------------
    // Default constructor: clear the dir_nodes, per_slave_nodes,
    // per_master_nodes; set num_dir_nodes, num_per_nodes to be zero;
    // set ID based on the above "no-nodal bc" setting.
    // --------------------------------------------------------------
    NodalBC( const int &nFunc );
    
    // --------------------------------------------------------------
    // The list of vtp files specifies the Dirichlet nodes. 
    // No periodical type BC nodes.
    // --------------------------------------------------------------
    NodalBC( const std::vector<std::string> &vtkfileList, const int &nFunc );

    // --------------------------------------------------------------
    // The list of vtk files specifies the Dirichlet nodes.
    // The list of sla files specifiesand periodical type BC nodes.
    // --------------------------------------------------------------
    NodalBC( const std::vector<std::string> &vtkfileList, const std::vector<std::string> &slafileList, const int &nFunc );

    // --------------------------------------------------------------
    // The list of vtp files specifies the rotated Dirichlet nodes. 
    // No periodical type BC nodes.
    // --------------------------------------------------------------
    NodalBC( const std::vector<std::string> &vtkfileList,
        const std::string &rotated_file, 
        const std::string &fixed_file,    
        const int &nFunc );

    // --------------------------------------------------------------
    // The list of vtp files specifies the rotated wall nodes. 
    // No periodical type BC nodes.
    // --------------------------------------------------------------
    NodalBC( const std::vector<std::string> &vtkfileList,
        const std::string &rotated_file,
        const std::string &rotated_wall_file, 
        const std::string &fixed_file,    
        const int &nFunc );

    // --------------------------------------------------------------
    // General NodalBC constructor. The specific implementation is 
    // determined by the private funcitons, indicated by the type var.
    // --------------------------------------------------------------
    NodalBC( const std::vector<std::string> &vtkfileList,
        const int &nFunc, const int &type );

    virtual ~NodalBC() = default;

    virtual unsigned int get_dir_nodes(const unsigned int &ii) const
    {return dir_nodes[ii];}

    virtual unsigned int get_per_slave_nodes(const unsigned int &ii) const
    {return per_slave_nodes[ii];}

    virtual unsigned int get_per_master_nodes(const unsigned int &ii) const
    {return per_master_nodes[ii];}

    virtual unsigned int get_num_dir_nodes() const {return num_dir_nodes;}

    virtual unsigned int get_num_per_nodes() const {return num_per_nodes;}

  private:
    std::vector<unsigned int> dir_nodes {};
    unsigned int num_dir_nodes {};

    std::vector<unsigned int> per_slave_nodes {}, per_master_nodes {};
    unsigned int num_per_nodes {};
    
    NodalBC() = delete;

    // --------------------------------------------------------------
    // BC_type_1 is a master-slave driven implementation. It is 
    // identical to the NodalBC(vtkfilelist, nFunc, master_idx),
    // with master_idx = [0,0,...,0].
    // --------------------------------------------------------------
    void BC_type_1( const std::vector<std::string> &vtkfileList,
        const int &nFunc );

    // --------------------------------------------------------------
    // BC_type_2 is a bc that sets Dirichlet on the first vtp surface
    // and set master-slave on the second vtp surface with master
    // node the node 0 in that vtp file.
    // This is used for the tensile test for arterial strips,
    // the first vtp file specifies the surface where the strip
    // is fixed; the second file specifies the surface where the
    // pulling force is applied.
    // --------------------------------------------------------------
    void BC_type_2( const std::vector<std::string> &vtkfileList,
        const int &nFunc );

    // --------------------------------------------------------------
    // BC_type_3 is a bc that sets Dirichlet on the given vtp file's
    // first node only. This is useful for Stokes/Navier-Stokes 
    // pressure boundary condition.
    // Note: the input vtifileList should have size 1, which means
    // only one vtp file is provided.
    // --------------------------------------------------------------
    void BC_type_3( const std::vector<std::string> &vtkfileList,
        const int &nFunc );
};

#endif
