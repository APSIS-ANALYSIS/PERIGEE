#ifndef EBC_PARTITION_OUTFLOW_HPP
#define EBC_PARTITION_OUTFLOW_HPP
// ==================================================================
// EBC_Partition_outflow.hpp
//
// Element boundary condition partition for outflow type boundary
// conditions.
// 
// If the partition owns the outlet face, we will record the basis 
// functions' face integral, and the outward normal vector.
//
// These information will be used for blood flow outlet BC setup.
// 
// Author: Ju Liu
// Date: Mar. 18 2019
// ==================================================================
#include "EBC_Partition.hpp"
#include "INodalBC.hpp"

class EBC_Partition_outflow : public EBC_Partition
{
  public:
    // The input ElemBC should be ElemBC_3D_outflow
    EBC_Partition_outflow( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const ElemBC * const &ebc,
        const std::vector<INodalBC *> &nbc_list );

    virtual ~EBC_Partition_outflow() = default;

    // write the data to hdf5 file in group /ebc/ebcid_xxx, 
    // xxx is the ebc_id
    virtual void write_hdf5( const std::string &FileName ) const
    { write_hdf5(FileName, "/ebc"); }

  protected:
    // Length is num_ebc x [ 0 if this part does not own this bc,
    // or ebc->get_num_node(ii) if this part owns this bc.] 
    std::vector< std::vector<double> > face_int_NA;

    // Length is num_ebc x [ 0 if this part does not own this bc,
    // or ebc -> get_num_node(ii) x 3 if this part owns this bc.]
    // 3 representing the nodes' x, y, z- velocity components.
    // both face_int_NA and LID_all_face_nodes are listed based on
    // the original surface vtp file node list, meaning
    // face_int_NA[][ii] and LID_all_face_nodes[][ii] are for the ii-th
    // nodes in the vtp file
    std::vector< std::vector<int> > LID_all_face_nodes;

    // unit normal vector to the cap surfaces, if the partition owns the bc.
    // size is num_ebc x [ 0, if does not own, or 3, if owns ].
    std::vector< Vector_3 > outvec; 
    
    // write the data to hdf5 file in group /group-name/ebcid_xxx, 
    // xxx is the ebc_id
    // We do not give users the access to this function out of the class
    virtual void write_hdf5( const std::string &FileName, 
        const std::string &GroupName ) const;

};

#endif
