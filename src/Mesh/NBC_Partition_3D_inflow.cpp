#include "NBC_Partition_3D_inflow.hpp"

NBC_Partition_3D_inflow::NBC_Partition_3D_inflow(
    const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const INodalBC * const &nbc ) : NBC_Partition_3D( part,
      mnindex, nbc )
{
  actarea  = nbc -> get_para_1();
  facearea = nbc -> get_para_6();

  outvec.clear();
  outvec.push_back( nbc->get_para_2(0) );
  outvec.push_back( nbc->get_para_2(1) );
  outvec.push_back( nbc->get_para_2(2) );

  num_out_bc_pts = nbc->get_para_3();

  centroid.resize(3);
  centroid[0] = nbc->get_para_4(0);
  centroid[1] = nbc->get_para_4(1);
  centroid[2] = nbc->get_para_4(2);

  outline_pts.resize( 3*num_out_bc_pts );
  for(int ii=0; ii<3*num_out_bc_pts; ++ii)
    outline_pts[ii] = nbc->get_para_5(ii);
}


NBC_Partition_3D_inflow:: ~NBC_Partition_3D_inflow()
{}


void NBC_Partition_3D_inflow::write_hdf5( const char * FileName ) const
{
  std::string filebname(FileName);
  std::string fName = SYS_T::gen_partfile_name( filebname, cpu_rank );
  
  hid_t file_id, group_id;

  file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  group_id = H5Gcreate(file_id, "/inflow", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5writer = new HDF5_Writer(file_id);

  if(LDN.size() > 0)
    h5writer->write_intVector(group_id, "LDN", LDN);

  h5writer->write_intVector(group_id, "Num_LD", Num_LD);

  h5writer->write_doubleVector( group_id, "Outward_normal_vector", outvec );

  h5writer->write_doubleScalar( group_id, "Inflow_active_area", actarea );
  
  h5writer->write_doubleScalar( group_id, "Inflow_full_area", facearea );

  h5writer->write_intScalar( group_id, "num_out_bc_pts", num_out_bc_pts);

  h5writer->write_doubleVector( group_id, "centroid", centroid);

  h5writer->write_doubleVector( group_id, "outline_pts", outline_pts);

  delete h5writer; H5Gclose(group_id); H5Fclose(file_id);
}

// EOF
