#ifndef AGLOBAL_MESH_INFO_FEM_3D_HPP
#define AGLOBAL_MESH_INFO_FEM_3D_HPP
// ==================================================================
// AGlobal_Mesh_Info_FEM_3D.hpp
//
// This is the global mesh info class used for a 3D FEM mesh.
//
// Date Created: Jan 21 2017
// ==================================================================
#include "IAGlobal_Mesh_Info.hpp"
#include "HDF5_Reader.hpp"

class AGlobal_Mesh_Info_FEM_3D final : public IAGlobal_Mesh_Info
{
  public:
    AGlobal_Mesh_Info_FEM_3D( const std::string &fileBaseName,
        const int &cpu_rank )
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
      elemType = h5r -> read_intScalar("Global_Mesh_Info", "elemType");

      H5Fclose( file_id );
    }

    ~AGlobal_Mesh_Info_FEM_3D() override = default;

    // Get the polynomial degree for the discretization method. For
    // unstructured mesh, the three function should return the same
    // value since one cannot differentiate the xyz direction in, e.g.,
    // tet mesh.
    int get_xdegree() const override {return xdegree;}
    int get_ydegree() const override {return ydegree;}
    int get_zdegree() const override {return zdegree;}

    // Get the total number of element of the whole mesh. 
    int get_nElem() const override {return nElem;}

    // Get the total number of nodes of the whole mesh.
    int get_nFunc() const override {return nFunc;}

    // Get the number of local basis functions for the element.
    // Note: this implicitly implies that we use the same type of element
    //       for the mesh.
    int get_nLocBas() const override {return nLocBas;}

    // Get the dimension of the problem.
    int get_probDim() const override {return probDim;}

    // Get an integer that indicate the element type.
    int get_elemType() const override {return elemType;}

    void print_info() const override
    {
      std::cout<<"AGlobal_Mesh_Info_FEM_3D:"<<std::endl;
      std::cout<<"degree: "<<xdegree<<'\t'<<ydegree<<'\t'<<zdegree<<'\n';
      std::cout<<"nElem: "<<nElem<<'\n';
      std::cout<<"nFunc: "<<nFunc<<'\n';
      std::cout<<"nLocBas: "<<nLocBas<<std::endl;
      std::cout<<"probDim: "<<probDim<<std::endl;
      std::cout<<"elemType: "<<get_elemType()<<std::endl;
    }

  private:
    int xdegree, ydegree, zdegree;
    int nElem, nFunc, nLocBas, probDim, elemType;

    AGlobal_Mesh_Info_FEM_3D() = delete;
};

#endif
