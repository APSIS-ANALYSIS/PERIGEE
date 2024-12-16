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

      const auto vdeg = h5r -> read_intVector("Global_Mesh_Info", "degree");

      xdegree = vdeg[0]; ydegree = vdeg[1]; zdegree = vdeg[2];

      nElem    = h5r -> read_intScalar("Global_Mesh_Info", "nElem");
      nFunc    = h5r -> read_intScalar("Global_Mesh_Info", "nFunc");
      nLocBas  = h5r -> read_intScalar("Global_Mesh_Info", "nLocBas");
      probDim  = h5r -> read_intScalar("Global_Mesh_Info", "probDim");
      elemType_str = h5r -> read_string("Global_Mesh_Info", "elemType");
      elemType = FE_T::to_FEType(elemType_str);

      H5Fclose( file_id );
    }

    ~AGlobal_Mesh_Info() = default;

    // --------------------------------------------------------------
    // Get the polynomial degree for the discretization method. For
    // unstructured mesh, the three function should return the same
    // value since one cannot differentiate the xyz direction in, e.g.,
    // tet mesh.
    // --------------------------------------------------------------
    int get_xdegree() const {return xdegree;}
    int get_ydegree() const {return ydegree;}
    int get_zdegree() const {return zdegree;}

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
      std::cout<<"degree: "<<xdegree<<'\t'<<ydegree<<'\t'<<zdegree<<'\n';
      std::cout<<"nElem: "<<nElem<<'\n';
      std::cout<<"nFunc: "<<nFunc<<'\n';
      std::cout<<"nLocBas: "<<nLocBas<<std::endl;
      std::cout<<"probDim: "<<probDim<<std::endl;
      std::cout<<"elemType: "<<elemType_str<<std::endl;
    }

  private:
    int xdegree, ydegree, zdegree;
    int nElem, nFunc, nLocBas, probDim;
    FEType elemType;
    std::string elemType_str;

    AGlobal_Mesh_Info() = delete;
};

#endif
