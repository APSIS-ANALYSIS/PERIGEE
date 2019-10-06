#include "EBC_Partition_2D.hpp"

EBC_Partition_2D::EBC_Partition_2D( const IPart * const &part,
            const Map_Node_Index * const &mnindex,
            const IElemBC * const &ebc )
: cpu_rank(part->get_cpu_rank())
{
  LFront_Elem.clear();
  LBack_Elem.clear();
  LLeft_Elem.clear();
  LRight_Elem.clear();
  LTop_Elem.clear();
  LBottom_Elem.clear();
  
  Num_LBCElem.clear();
  Num_LBCElem.resize(6); // In general, there are 6 faces
  
  Num_LBCElem[4] = 0; // top 
  Num_LBCElem[5] = 0; // bottom

  unsigned int lfront_elem = 0;
  int temp_elem, elemlocindex;
  for(unsigned int jj=0; jj<ebc->get_num_front_elem(); ++jj)
  {
    temp_elem = ebc->get_front_elem(jj);
    elemlocindex = part->get_elemLocIndex(temp_elem);
    if(elemlocindex != -1)
    {
      LFront_Elem.push_back(elemlocindex);
      lfront_elem += 1;
    }
  }
  Num_LBCElem[0] = lfront_elem; // front
  
  
  unsigned int lback_elem = 0;
  for(unsigned int jj=0; jj<ebc->get_num_back_elem(); ++jj)
  {
    temp_elem = ebc->get_back_elem(jj);
    elemlocindex = part->get_elemLocIndex(temp_elem);
    if(elemlocindex != -1)
    {
      LBack_Elem.push_back(elemlocindex);
      lback_elem += 1;
    }
  } 
  Num_LBCElem[1] = lback_elem; // back

  unsigned int lleft_elem = 0;
  for(unsigned int jj=0; jj<ebc->get_num_left_elem(); ++jj)
  {
    temp_elem = ebc->get_left_elem(jj);
    elemlocindex = part->get_elemLocIndex(temp_elem);
    if(elemlocindex != -1)
    {
      LLeft_Elem.push_back(elemlocindex);
      lleft_elem += 1;
    }
  }
  Num_LBCElem[2] = lleft_elem; // left

  unsigned int lrig_elem = 0;
  for(unsigned int jj=0; jj<ebc->get_num_right_elem(); ++jj)
  {
    temp_elem = ebc->get_right_elem(jj);
    elemlocindex = part->get_elemLocIndex(temp_elem);
    if(elemlocindex != -1)
    {
      LRight_Elem.push_back(elemlocindex);
      lrig_elem += 1;
    }
  }
  Num_LBCElem[3] = lrig_elem; // right


  VEC_T::shrink2fit(LFront_Elem);
  VEC_T::shrink2fit(LBack_Elem);
  VEC_T::shrink2fit(LLeft_Elem);
  VEC_T::shrink2fit(LRight_Elem);
}



EBC_Partition_2D::~EBC_Partition_2D()
{
  VEC_T::clean(LFront_Elem);
  VEC_T::clean(LBack_Elem);
  VEC_T::clean(LLeft_Elem);
  VEC_T::clean(LRight_Elem);
  VEC_T::clean(LTop_Elem);
  VEC_T::clean(LBottom_Elem);
  VEC_T::clean(Num_LBCElem);
}



void EBC_Partition_2D::write_hdf5(const char * FileName) const
{
  std::string fName(FileName);
  fName.append("_p");

  if( cpu_rank / 10 == 0 )
    fName.append("0000");
  else if( cpu_rank / 100 == 0 )
    fName.append("000");
  else if( cpu_rank / 1000 == 0 )
    fName.append("00");
  else if( cpu_rank / 10000 == 0 )
    fName.append("0");

  std::stringstream sstrm;
  sstrm<<cpu_rank;
  fName.append(sstrm.str());

  fName.append(".h5");

  hid_t file_id, group_id;

  file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  group_id = H5Gcreate(file_id, "/ebc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5writer = new HDF5_Writer(file_id);

  h5writer->write_intVector(group_id, "Num_LBCElem", Num_LBCElem);

  if(LLeft_Elem.size() > 0)
    h5writer->write_intVector(group_id, "LLeft_Elem", LLeft_Elem);

  if(LRight_Elem.size() > 0)
    h5writer->write_intVector(group_id, "LRight_Elem", LRight_Elem);

  if(LFront_Elem.size() > 0)
    h5writer->write_intVector(group_id, "LFront_Elem", LFront_Elem);

  if(LBack_Elem.size() > 0)
    h5writer->write_intVector(group_id, "LBack_Elem", LBack_Elem);



  delete h5writer;
  H5Gclose(group_id);
  H5Fclose(file_id);
}



void EBC_Partition_2D::print_info() const
{
  std::cout<<"=========================================== \n";
  std::cout<<"EBC_Partition_2D : \n";
  std::cout<<"--- Num_LBCElem : ";
  VEC_T::print(Num_LBCElem);
  std::cout<<"\n--- LFront_Elem : \n";
  VEC_T::print(LFront_Elem);
  std::cout<<"\n--- LBack_Elem : \n";
  VEC_T::print(LBack_Elem);
  std::cout<<"\n--- LLeft_Elem : \n";
  VEC_T::print(LLeft_Elem);
  std::cout<<"\n--- LRight_Elem : \n";
  VEC_T::print(LRight_Elem);
  std::cout<<"\n--- LTop_Elem : \n";
  VEC_T::print(LTop_Elem);
  std::cout<<"\n--- LBottom_Elem : \n";
  VEC_T::print(LBottom_Elem);
  std::cout<<"=========================================== \n";
}

// EOF
