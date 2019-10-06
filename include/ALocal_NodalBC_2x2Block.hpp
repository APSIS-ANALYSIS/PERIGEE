#ifndef ALOCAL_NODALBC_2X2BLOCK_HPP
#define ALOCAL_NODALBC_2X2BLOCK_HPP

class ALocal_NodalBC_2x2Block
{
  public:
    ALocal_NodalBC_2x2Block(){};

    virtual ~ALocal_NodalBC_2x2Block(){};

    virtual void print_info() const = 0;

    virtual int get_LID_0( const int &node ) const {return -1;}

    virtual int get_LID_1( const int &dof_idx, const int &node ) const
    {return -1;}

    virtual int get_LDN_0( const int &node ) const {return -1;}

    virtual int get_LDN_1( const int &dof_idx, const int &node ) const 
    {return -1;}

    virtual int get_Num_LD_0() const {return -1;}

    virtual int get_Num_LD_1( const int &dof_idx ) const {return -1;}
};


#endif
