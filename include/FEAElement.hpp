#ifndef FEAELEMENT_HPP
#define FEAELEMENT_HPP
// ============================================================================
// FEAElement.hpp
// Object: This is the interface for Finite Element classes.
// It records:
// 1. Element basic indices;
// 2. Element nodes through the IEN array;
// 3. reference element position / hx hy hz; 
//
// Author: Ju Liu
// Date created: Nov. 6 2013
// ============================================================================
#include "IQuadPts.hpp"
#include "FEANode.hpp"

enum class ElementType
{
	Tet4,
	Tet10,
	Hex8,
	Hex27,
};

class FEAElement
{
  public:
    // Constructor
    FEAElement() = default;
    
    // Destructor
    virtual ~FEAElement() = default;

    // Return this element's dimension
    virtual int get_elemDim() const = 0;

    // Return this element's Type, which defines the type of different 
    // elements defined on this single element domain.
    virtual int get_Type() const
    {SYS_T::commPrint("Warning: get_Type is not implemented. \n"); return -1;}

    // Return the element's type name
    virtual std::string get_TypeName() const
    {SYS_T::commPrint("Warning: get_TypeName is not implemented. \n"); return "undetermined";}

    // Return the number of nodes for each type elements
    virtual int get_nLocBas() const = 0;
    
    // Return the quadrature info
    virtual int get_numQuapts() const = 0;

    // ------------------------------------------------------------------------
    // print function
    // ------------------------------------------------------------------------
    virtual void print_info() const 
    {SYS_T::commPrint("Warning: print is not implemented. \n");}

    // Return the memory usage of this class in bytes
    virtual double get_memory_usage() const = 0; 

    // ------------------------------------------------------------------------
    // Calculate the element size
    // ------------------------------------------------------------------------
    virtual double get_h( const double * const &ctrl_x, const double * const &ctrl_y,
        const double * const &ctrl_z ) const
    {SYS_T::commPrint("Warning: get_h is not implemented. \n"); return 0.0;}

    virtual double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y ) const
    {SYS_T::commPrint("Warning: get_h is not implemented. \n"); return 0.0;}

    virtual double get_h( const std::array<std::vector<double>,3> pts ) const
    {return get_h( &pts[0][0], &pts[1][0], &pts[2][0] );}

    virtual double get_h( const std::array<std::vector<double>,2> pts ) const
    {return get_h( &pts[0][0], &pts[1][0] );}

    // ------------------------------------------------------------------------
    // Build Basis function quadrature info
    // ------------------------------------------------------------------------
    // Build 3D basis -- FEM
    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x, const double * const &ctrl_y,
        const double * const &ctrl_z )
    {SYS_T::commPrint("Warning: buildBasis() is not implemented. \n");}

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const std::array<std::vector<double>,3> &pts )
    {buildBasis(quad_rule, &pts[0][0], &pts[1][0], &pts[2][0]);}

    // Build 2D basis -- FEM
    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x, const double * const &ctrl_y )
    {SYS_T::commPrint("Warning: buildBasis() is not implemented. \n");}

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const std::array<std::vector<double>,2> &pts )
    {buildBasis(quad_rule, &pts[0][0], &pts[1][0]);}

    // ------------------------------------------------------------------------
    // Get functions : Obtain the value of basis functions and their derivatives.
    // ------------------------------------------------------------------------
    // Return basis function value using a dynamic double array.
    // Users are responsible for allocating proper memory for basis, and delete 
    // the pointer after use.
    // ------------------------------------------------------------------------
    virtual void get_R( const int &quaindex, double * const &basis ) const
    {SYS_T::commPrint("Warning: get_R is not implemented. \n");} 

    virtual std::vector<double> get_R( const int &quaindex ) const
    {SYS_T::commPrint("Warning: get_R is not implemented. \n"); return {};} 

    // ------------------------------------------------------------------------    
    // Return basis function and all 1st order derivatives value using a dynamic
    // array.
    // Users are responsible for allocating proper memory for basis, and delete
    // the pointer after use.
    // ------------------------------------------------------------------------    
    // 3D case
    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y, double * const &basis_z ) const 
    {SYS_T::commPrint("Warning: get_gradR is not implemented. \n");} 

    // 2D case
    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y ) const 
    {SYS_T::commPrint("Warning: get_gradR is not implemented. \n");} 

    // 2D case:
    virtual void get_R_gradR( const int &quaindex, double * const &basis, 
        double * const &basis_x, double * const &basis_y ) const 
    {SYS_T::commPrint("Warning: get_R_gradR is not implemented. \n");} 

    // 3D case:
    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y, double * const &basis_z ) const 
    {SYS_T::commPrint("Warning: get_R_gradR is not implemented. \n");}

    virtual std::vector<double> get_dR_dx( const int &quaindex ) const
    {SYS_T::commPrint("Warning: get_dR_dx is not implemented. \n");return {};}

    virtual std::vector<double> get_dR_dy( const int &quaindex ) const
    {SYS_T::commPrint("Warning: get_dR_dy is not implemented. \n");return {};}

    virtual std::vector<double> get_dR_dz( const int &quaindex ) const
    {SYS_T::commPrint("Warning: get_dR_dz is not implemented. \n");return {};}

    // ------------------------------------------------------------------------    
    // R, gradR, and Laplacian R
    // ------------------------------------------------------------------------    
    virtual void get_3D_R_gradR_LaplacianR( const int &quaindex,
        double * const &basis, double * const &basis_x, double * const &basis_y,
        double * const &basis_z, double * const &basis_xx, double * const &basis_yy, 
        double * const &basis_zz ) const 
    {SYS_T::commPrint("Warning: get_3DLaplacianR is not implemented. \n");}

    virtual void get_2D_R_gradR_LaplacianR( const int &quaindex,
        double * const &basis, double * const &basis_x, double * const &basis_y,
        double * const &basis_xx, double * const &basis_yy ) const 
    {SYS_T::commPrint("Warning: get_2DLaplacianR is not implemented. \n");}

    // ------------------------------------------------------------------------    
    // R, gradR, and grad gradR
    // ------------------------------------------------------------------------    
    virtual void get_2D_R_dR_d2R( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y, double * const &basis_xx, 
        double * const &basis_yy, double * const &basis_xy ) const
    {SYS_T::commPrint("Warning: get_2D_R_dR_d2R is not implemented. \n");}

    virtual void get_3D_R_dR_d2R( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y, double * const &basis_z,
        double * const &basis_xx, double * const &basis_yy, double * const &basis_zz,
        double * const &basis_xy, double * const &basis_xz, double * const &basis_yz ) 
      const {SYS_T::commPrint("Warning: get_3D_R_dR_d2R is not implemented. \n");}

    virtual std::vector<double> get_d2R_dxx( const int &quaindex ) const
    {SYS_T::commPrint("Warning: get_d2R_dxx is not implemented. \n");return {};}

    virtual std::vector<double> get_d2R_dyy( const int &quaindex ) const
    {SYS_T::commPrint("Warning: get_d2R_dyy is not implemented. \n");return {};}

    virtual std::vector<double> get_d2R_dzz( const int &quaindex ) const
    {SYS_T::commPrint("Warning: get_d2R_dzz is not implemented. \n");return {};}

    virtual std::vector<double> get_d2R_dxy( const int &quaindex ) const
    {SYS_T::commPrint("Warning: get_d2R_dxy is not implemented. \n");return {};}

    virtual std::vector<double> get_d2R_dxz( const int &quaindex ) const
    {SYS_T::commPrint("Warning: get_d2R_dxz is not implemented. \n");return {};}

    virtual std::vector<double> get_d2R_dyz( const int &quaindex ) const
    {SYS_T::commPrint("Warning: get_d2R_dyz is not implemented. \n");return {};}

    // ------------------------------------------------------------------------    
    // Return the Jacobian determinant
    // ------------------------------------------------------------------------    
    virtual double get_detJac(const int &quaindex) const = 0;

    // ------------------------------------------------------------------------    
    // Return the Jacobian matrix in rows, i.e. dx_dxi, dx_deta, dx_dzeta,
    // dy_dxi, ... dz_deta, dz_dzeta. The size of teh return vector is 9 for 3d,
    // and 4 for 2d, for exampel. 
    // ------------------------------------------------------------------------    
    virtual void get_Jacobian(const int &quaindex, double * const &jac_value) const 
    {SYS_T::commPrint("Warning: get_Jacobian is not implemented. \n");} 

    virtual std::array<double,9> get_Jacobian( const int &quaindex ) const
    {
      SYS_T::commPrint("Warning: get_Jacobian is not implemented. \n");
      return {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    }

    // ------------------------------------------------------------------------    
    // Return the inversion of the Jacobian matrix, the 9(4) components are
    // given. The output array dxi_dx is a 9 or 4 component array for 3D / 2D
    // case.
    // ------------------------------------------------------------------------    
    virtual void get_invJacobian(const int &quaindex, double * const &dxi_dx) const 
    {SYS_T::commPrint("Warning: get_invJacobian is not implemented. \n");} 

    virtual std::array<double,9> get_invJacobian( const int &quaindex ) const
    {
      SYS_T::commPrint("Warning: get_Jacobian is not implemented. \n");
      return {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    }

    // ------------------------------------------------------------------------    
    // Return the global-to-lamina rotation matrix for transforming between the
    // global and local lamina coordinates. This routine is used in membrane
    // and shell elements.
    // Reference: TJRH Linear finite element book page 386.
    // ------------------------------------------------------------------------    
    virtual Tensor2_3D get_rotationMatrix( const int &quaindex ) const 
    {
      SYS_T::commPrint("Warning: get_rotationMatrix is not implemented. \n"); 
      return Tensor2_3D();
    }

    // ------------------------------------------------------------------------
    // Unit outward normal vector and boundary area measure.
    // In a oriented surface, the outward normal direction and the area
    // measure can be calculated using the cross product of vector x_,xi and
    // vector x_,eta. Assuing that the nodal points are properly ordered, this
    // product will point in the outward direction with respect to the volume
    // element. To ensure the ordering, some special procedure may need to be
    // performed in the preprocessor to check and order the nodal indices.
    // See, ElemBC_3D::resetSurIEN_outwardnormal() function for example
    // for the ordering of a triangle element that belongs to a tetrahedron
    // surface.
    // This function is called in FEAElement_Triangle3_3D_der0.
    // ------------------------------------------------------------------------
    virtual Vector_3 get_2d_normal_out( const int &quaindex, double &area ) const
    {
      SYS_T::commPrint("Warning: get_2d_normal_out is not implemented. \n");
      return Vector_3();
    }

    // ------------------------------------------------------------------------
    // Unit outward normal vector and boundary area measure.
    // In a line element, the tangential vector is well defined. Given
    // an interior point's coordinate, the outward normal vector is defined.
    // intpt_x/y/z gives the interior point's coordinate.
    // quaindex gives the index of the quadrature point, from which one can
    // calculate the "root" of the tangential vector.
    // This function is, for example, called in FEAElement_Line2_3D_der0.
    // ------------------------------------------------------------------------
    virtual Vector_3 get_normal_out( const int &quaindex,
        const std::vector< Vector_3> &sur_pt, const Vector_3 &int_pt, 
        double &length ) const
    {
      SYS_T::commPrint("Warning: get_normal_out is not implemented. \n");
      return Vector_3();
    }

    virtual Vector_3 get_normal_out( const int &quaindex,
        const Vector_3 &sur_pt, const Vector_3 &int_pt, 
        double &length ) const
    {
      SYS_T::commPrint("Warning: get_normal_out is not implemented. \n");
      return Vector_3();
    }

    // ------------------------------------------------------------------------
    // Build the volume element with a face id and the quad_rule on surface element.
    // ------------------------------------------------------------------------
    virtual void buildBasis( const int &face_id,
        const IQuadPts * const &quad_rule_s,
        const double * const &ctrl_x, 
        const double * const &ctrl_y,
        const double * const &ctrl_z )
    {SYS_T::commPrint("Warning: buildBasis is not implemented. \n");}

    // ------------------------------------------------------------------------
    // dx_dr in parent domain
    // ------------------------------------------------------------------------
    virtual Vector_3 get_dx_dr( const int &quaindex,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const 
    {
      SYS_T::commPrint("Warning: get_dx_dr is not implemented. \n");
      return Vector_3();
    }

    // ------------------------------------------------------------------------
    // dx_ds in parent domain
    // ------------------------------------------------------------------------
    virtual Vector_3 get_dx_ds( const int &quaindex,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const 
    {
      SYS_T::commPrint("Warning: get_dx_ds is not implemented. \n");
      return Vector_3();
    }

    // ------------------------------------------------------------------------
    // d2x_drr in parent domain
    // ------------------------------------------------------------------------
    virtual Vector_3 get_d2x_drr( const int &quaindex,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const 
    {
      SYS_T::commPrint("Warning: get_d2x_drr is not implemented. \n");
      return Vector_3();
    }

    // ------------------------------------------------------------------------
    // d2x_dss in parent domain
    // ------------------------------------------------------------------------
    virtual Vector_3 get_d2x_dss( const int &quaindex,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const 
    {
      SYS_T::commPrint("Warning: get_d2x_dss is not implemented. \n");
      return Vector_3();
    }

    // ------------------------------------------------------------------------
    // d2x_drs in parent domain
    // ------------------------------------------------------------------------
    virtual Vector_3 get_d2x_drs( const int &quaindex,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const 
    {
      SYS_T::commPrint("Warning: get_d2x_drs is not implemented. \n");
      return Vector_3();
    }

    virtual std::array<std::vector<double>, 3> get_face_ctrlPts( const int &face_id,
        const double * const &volctrl_x,
        const double * const &volctrl_y,
        const double * const &volctrl_z )
    {
      SYS_T::commPrint("Warning: get_face_ctrlPts is not implemented. \n");
    
      return {{ }};
    }

};

#endif
