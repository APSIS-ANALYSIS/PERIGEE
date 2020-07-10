#ifndef EBC_PARTITION_VTP_OUTFLOW_HPP
#define EBC_PARTITION_VTP_OUTFLOW_HPP
// ==================================================================
// EBC_Partition_vtp_outflow.hpp
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
#include "EBC_Partition_vtp.hpp"
#include "INodalBC.hpp"

class EBC_Partition_vtp_outflow : public EBC_Partition_vtp
{
  public:
    // The input ElemBC should be ElemBC_3D_tet4_outflow
    EBC_Partition_vtp_outflow( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const ElemBC * const &ebc,
        const std::vector<INodalBC *> &nbc_list );

    virtual ~EBC_Partition_vtp_outflow();

    // write the data to hdf5 file in group /ebc/ebcid_xxx, 
    // xxx is the ebc_id
    virtual void write_hdf5( const char * FileName ) const;

    // write the data to hdf5 file in group /group-name/ebcid_xxx, 
    // xxx is the ebc_id
    virtual void write_hdf5( const char * FileName, const char * GroupName ) const;

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

    // unit normal vector to the faces, if the partition owns the bc.
    // size is num_ebc x [ 0, if does not own, or 3, if nows ].
    std::vector< std::vector<double> > outvec; 
};

#endif
