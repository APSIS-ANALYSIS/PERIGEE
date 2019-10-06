#include "Analysis_Tool_Test.hpp"

void AnalysisTool_EXT_check(const std::string &part_file, const int &rank)
{
  HDF5_PartReader * h5reader = new HDF5_PartReader(part_file, rank);
 
  IAExtractor * fExt = new AExtractor_3D_NURBS_xyz(part_file, rank);

  int nelem;
  std::vector<int> elem_loc;

  h5reader->get_LE(elem_loc, nelem);

  for(int ee=0; ee<nelem; ++ee)
  {
    std::vector<double> ext_x, ext_y, ext_z;
    //h5reader->get_EXT_x(ee, ext_x);
    //h5reader->get_EXT_y(ee, ext_y);
    //h5reader->get_EXT_z(ee, ext_z);

    std::vector<double> ext_xx, ext_yy, ext_zz;
    //fExt->get_EXT_x(ee, ext_xx);
    //fExt->get_EXT_y(ee, ext_yy);
    //fExt->get_EXT_z(ee, ext_zz);
    
    assert(ext_x.size() == ext_xx.size());
    assert(ext_y.size() == ext_yy.size());
    assert(ext_z.size() == ext_zz.size());

    for(unsigned int ii=0; ii<ext_x.size(); ++ii)
      assert(ext_x[ii] == ext_xx[ii]);

    for(unsigned int ii=0; ii<ext_y.size(); ++ii)
      assert(ext_y[ii] == ext_yy[ii]);
    
    for(unsigned int ii=0; ii<ext_z.size(); ++ii)
      assert(ext_z[ii] == ext_zz[ii]);
  }

  std::cout<<"IAExtractor : Extractor_3D_NURBS_xyz is correct!\n";
  delete fExt;
  delete h5reader;
}

void AnalysisTool_LIEN_check(const HDF5_PartReader * const &h5reader)
{
  ALocal_IEN * alien = new ALocal_IEN(h5reader);
  
  int nlocalele;
  std::vector<int> temp;
  h5reader->get_LE(temp, nlocalele);
  temp.clear();

  std::vector<int> temp_1, temp_2;

  for(int ee=0; ee<nlocalele; ++ee)
  {
    h5reader->get_LIEN(ee, temp_1);
    alien->get_LIEN_e(ee, temp_2);

    assert(temp_1.size() == temp_2.size());
  
    for(unsigned int ii=0; ii<temp_1.size(); ++ii)
      assert(temp_1[ii] == temp_2[ii]);
  }

  std::cout<<"LIEN is read correctly.\n";
  delete alien;
}

