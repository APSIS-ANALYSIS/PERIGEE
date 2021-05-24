# A List of Changes

* PDNSolution.hpp

  - GetLocalArray( std::vector<double> ) function is changed to std::vector<double> GetLocalArray(). 79c7991cf5b9a0262adc5aa4bf0eb8c2030240c4

* ALocal_IEN.hpp
  - get_LIEN_e functions are renamed as get_LIEN. 9b6bcccf4ddb72a85fbb04a7316644912a0cd467
  - get_LIEN( const int &elem, std::vector<int> &elem_ien ) is changed to std::vector<int> geT_LIEN( const int &elem ).
