#ifndef AGLOBAL_MESH_INFO_HPP
#define AGLOBAL_MESH_INFO_HPP
// ==================================================================
// AGlobal_Mesh_Info.hpp
//
// This is the global mesh info class used for a 3D FEM mesh.
//
// Date Created: Jan 21 2017
// ==================================================================
#include "HDF5_Reader.hpp"
#include "FEType.hpp"

class AGlobal_Mesh_Info
{
  public:
    AGlobal_Mesh_Info( const std::string &fileBaseName, const int &cpu_rank )
    {
      const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

      hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

      std::unique_ptr<HDF5_Reader> h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

      nElem    = h5r -> read_intScalar("Global_Mesh_Info", "nElem");
      nFunc    = h5r -> read_intScalar("Global_Mesh_Info", "nFunc");
      nLocBas  = h5r -> read_intScalar("Global_Mesh_Info", "nLocBas");
      probDim  = h5r -> read_intScalar("Global_Mesh_Info", "probDim");
      elemType = FE_T::to_FEType(h5r -> read_string("Global_Mesh_Info", "elemType"));

      H5Fclose( file_id );
    }

    AGlobal_Mesh_Info( const HDF5_Reader * const &h5r )
    {
      nElem    = h5r -> read_intScalar("Global_Mesh_Info", "nElem");
      nFunc    = h5r -> read_intScalar("Global_Mesh_Info", "nFunc");
      nLocBas  = h5r -> read_intScalar("Global_Mesh_Info", "nLocBas");
      probDim  = h5r -> read_intScalar("Global_Mesh_Info", "probDim");
      elemType = FE_T::to_FEType(h5r -> read_string("Global_Mesh_Info", "elemType"));
    }

    ~AGlobal_Mesh_Info() = default;

    // --------------------------------------------------------------
    // Get the total number of element of the whole mesh. 
    // --------------------------------------------------------------
    int get_nElem() const {return nElem;}

    // --------------------------------------------------------------
    // Get the total number of nodes of the whole mesh.
    // --------------------------------------------------------------
    int get_nFunc() const {return nFunc;}

    // --------------------------------------------------------------
    // Get the number of local basis functions for the element.
    // Note: this implicitly implies that we use the same type of element
    //       for the mesh.
    // --------------------------------------------------------------
    int get_nLocBas() const {return nLocBas;}

    // --------------------------------------------------------------
    // Get the dimension of the problem.
    // --------------------------------------------------------------
    int get_probDim() const {return probDim;}

    // --------------------------------------------------------------
    // Get an integer that indicate the element type.
    // --------------------------------------------------------------
    FEType get_elemType() const {return elemType;}

    void print_info() const
    {
      std::cout<<"AGlobal_Mesh_Info:"<<std::endl;
      std::cout<<"nElem: "<<nElem<<'\n';
      std::cout<<"nFunc: "<<nFunc<<'\n';
      std::cout<<"nLocBas: "<<nLocBas<<std::endl;
      std::cout<<"probDim: "<<probDim<<std::endl;
      std::cout<<"elemType: "<<FE_T::to_string(elemType)<<std::endl;
    }

  private:
    int nElem, nFunc, nLocBas, probDim;
    FEType elemType;

    AGlobal_Mesh_Info() = delete;
};

#endif
