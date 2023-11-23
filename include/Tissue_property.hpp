#ifndef TISSUE_PROPERTY_HPP
#define TISSUE_PROPERTY_HPP
// ============================================================================
// Tissue_property.hpp
//
// This class stores the direction basis vectors at solid nodes. 
// Direction basis includes 3 Vector_3 for each node.
//
// Constructed from inputs of the file basename and cpu_rank, and read using 
// the HDF5_Reader tool.
//
// Date: Apr. 13th 2023
// Author: Qingshuang Lu
// ============================================================================
#include "HDF5_Reader.hpp"

class Tissue_property
{
  public:
    Tissue_property( const std::string &fileBaseName, const int &rank );

    virtual ~Tissue_property() = default;

    // ------------------------------------------------------------------------
    // Input: nn the index of local_to_global, including local nodes and ghost nodes. 
    //        The maxmum of nn is nlocghonode-1.
    // Output: the Vector_3 of the corresponding direction basis vector
    // ------------------------------------------------------------------------
    virtual Vector_3 get_basis_r( const int &nn ) const;

    virtual Vector_3 get_basis_l( const int &nn ) const;

    virtual Vector_3 get_basis_c( const int &nn ) const;

    // ------------------------------------------------------------------------
    // ! Print the info for this class.
    // ------------------------------------------------------------------------
    virtual void print_info() const;

  private:
    // ------------------------------------------------------------------------
    // The rank or id of the subdomain
    // ------------------------------------------------------------------------
    const int cpu_rank;

    // ------------------------------------------------------------------------
    // The number of local solid node
    // ------------------------------------------------------------------------
    int nlocghonode_s;

    // ------------------------------------------------------------------------
    // The mapping from the total local nodes to the local solid nodes.
    // ------------------------------------------------------------------------
    std::vector<int> node_locgho_solid;

    // ------------------------------------------------------------------------
    // The direction basis vectors. Letters r, l, and c denote radial, longitudinal,
    // and circumferential, respectively. 
    // ------------------------------------------------------------------------
    std::vector< Vector_3 > basis_r, basis_l, basis_c;
};

#endif
