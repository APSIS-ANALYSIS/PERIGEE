#include "EBC_Partition_weak.hpp"

EBC_Partition_weak::EBC_Partition_weak(const IPart * const &part,
    const Map_Node_Index * const &mnindex, const ElemBC * const &ebc)
: EBC_Partition(part, mnindex, ebc),
weak_bc_type {ebc->get_weak_bc_type()}, C_bI {ebc->get_C_bI()}
{
  if(weak_bc_type == 0)
    ;   // do nothing
  else if(weak_bc_type == 1 || weak_bc_type == 2)
  {

    vol_ele_id.resize(num_ebc);
    ele_face_id.resize(num_ebc);
    for(int ii{0}; ii < num_ebc; ++ii)
    {   
        vol_ele_id.resize(get_num_local_cell(ii));
        ele_face_id.resize(get_num_local_cell(ii));
        for(int ee{0}; ee < get_num_local_cell(ii); ++ee)
        {
          vol_ele_id[ii][ee] = ebc->get_global_cell(ii, ee);
          ele_face_id[ii][ee] = ebc->get_face_id(ii, ee);
        }
    }
  }
  else
    SYS_T::print_fatal("Error: EBC_Partition_weak, unknown weak bc type.\n");

  if(weak_bc_type == 2)
    ;   // Rotation matrix, unimplemented
}

EBC_Partition_weak::~EBC_Partition_weak()
{
  vol_ele_id.clear();
  ele_face_id.clear();
}