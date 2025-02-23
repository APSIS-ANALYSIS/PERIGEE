#ifndef PLOCASSEM_FSI_MESH_LAPLACIAN_HPP
#define PLOCASSEM_FSI_MESH_LAPLACIAN_HPP
// ==================================================================
// PLocAssem_FSI_Mesh_Laplacian.hpp
//
// This is the local assembly routine for the mesh motion in FSI
// problem using the harmonic extension of the solid displacement.
//
// Date Created: July 31 2017
// Author: Ju Liu
// ==================================================================
#include "IPLocAssem.hpp"

class PLocAssem_FSI_Mesh_Laplacian : public IPLocAssem
{
  public:
    PLocAssem_FSI_Mesh_Laplacian( const FEType &in_type, const int &in_nqp_v, const int &in_nqp_s );

    virtual ~PLocAssem_FSI_Mesh_Laplacian();

    virtual int get_dof() const {return 3;}

    virtual int get_dof_mat() const {return 3;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual int get_snLocBas() const {return snLocBas;}

    virtual int get_num_ebc_fun() const {return num_ebc_fun;}

    virtual void Zero_Tangent_Residual();

    virtual void Zero_Residual();

    virtual void Assem_Estimate();

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    virtual void Assem_Mass_Residual(
        const double * const &vec_b,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {SYS_T::print_fatal("Error: PLocAssem_FSI_Mesh_Laplacian::Assem_Mass_Residual is not implemented. \n");}

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )
    {SYS_T::print_fatal("Error: PLocAssem_FSI_Mesh_Laplacian::Assem_Residual_EBC is not implemented. \n");}

  private:
    const FEType elemType;

    const int nqpv, nqps;

    const std::unique_ptr<FEAElement> elementv, elements;

    const std::unique_ptr<IQuadPts> quadv, quads;

    const int num_ebc_fun, vec_size;

    void print_info() const;

    Vector_3 get_f( const double &x, const double &y, const double &z,
        const double &t ) const
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    void get_ebc_g1( const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }

    typedef void ( PLocAssem_FSI_Mesh_Laplacian::*locassem_fsi_mesh_lap_funs )( const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const;

    locassem_fsi_mesh_lap_funs * flist;

    void get_ebc_fun( const int &ebc_id,
        const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      return ((*this).*(flist[ebc_id]))(x,y,z,t,nx,ny,nz,gx,gy,gz);
    }
};

#endif
