#ifndef PDNSOLUTION_WALL_DISP_HPP
#define PDNSOLUTION_WALL_DISP_HPP
// ==================================================================
// PDNSolution_Wall_Disp.hpp
//
// This is the solution vector for CMM solver -- wall displacement.
//
// The solution has 3 dofs, that is the three componenet of displacement.
// ==================================================================
#include "Math_Tools.hpp"
#include "PDNSolution.hpp"
#include "FEANode.hpp"
#include "ALocal_EBC.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"

class PDNSolution_Wall_Disp : public PDNSolution
{
  public:
    PDNSolution_Wall_Disp( const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const int &type, const bool &isprint = true );

    // ==== WOMERSLEY CHANGES BEGIN ====
    PDNSolution_Wall_Disp( const IAGlobal_Mesh_Info * const &agmi_ptr,
        const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const ALocal_EBC * const &ebc_wall_part,
        const double &rho,
        const int &type, const bool &isprint = true );
    // ==== WOMERSLEY CHANGES END ====

    ~PDNSolution_Wall_Disp();

    // case 0 : generate a full zero solution vector
    void Init_zero( const APart_Node * const &pNode_ptr );

    // ==== WOMERSLEY CHANGES BEGIN ====

    // case 1 : generate Womersley wall disp
    void Init_womersley( const IAGlobal_Mesh_Info * const &agmi_ptr,
        const APart_Node * const &pNode_ptr,
        const ALocal_EBC * const &ebc_wall_part,
        const double &rho );

    // case 2 : generate Womersley dot wall disp
    void Init_womersley_dot( const IAGlobal_Mesh_Info * const &agmi_ptr,
        const APart_Node * const &pNode_ptr,
        const ALocal_EBC * const &ebc_wall_part,
        const double &rho );

    void ReadNodeMapping( const char * const &node_mapping_file,
        const char * const &mapping_type, const int &node_size,
        int * const &nodemap );

    // ==== WOMERSLEY CHANGES END ====

  private:
    const bool is_print;
};

#endif
