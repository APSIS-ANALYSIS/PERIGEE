# A List of Changes

* PDNSolution.hpp

  - GetLocalArray( std::vector<double> ) function is changed to std::vector<double> GetLocalArray(). 79c7991cf5b9a0262adc5aa4bf0eb8c2030240c4

* ALocal_IEN.hpp
  - get_LIEN_e functions are renamed as get_LIEN. 9b6bcccf4ddb72a85fbb04a7316644912a0cd467
  - get_LIEN( const int &elem, std::vector<int> &elem_ien ) is changed to std::vector<int> geT_LIEN( const int &elem ). d4957932e3d21896eef89da991af2c536c6689df

* FEAElement.hpp
  - get_2d_normal_out function is changed to return the normal vector as a Vector_3 object directly. 057cc3c93d88359ae35e98eff6688809423d6743
