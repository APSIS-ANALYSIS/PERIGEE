#include "QuadPts_Gauss_Tet_face.hpp"

QuadPts_Gauss_Tet_face::QuadPts_Gauss_Tet_face( const int &in_num_pts, const int &in_face_id )
: num_pts( in_num_pts ), face_id( in_face_id ), QP_triangle( QuadPts_Gauss_Triangle( in_num_pts ))
{
  SYS_T::print_fatal_if( face_id < 0 || face_id > 3,
    "Error: QuadPts_Gauss_Tet_face, wrong face id input.\n");

  qp.assign( 4 * num_pts, 0.0 );

  switch(face_id)
  {
    case 0: // u = 0 : node1 = node0', node2 = node1', node3 = node2'
      for(unsigned int ii{0}; ii < num_pts; ++ii)
      {
        qp[4*ii + 0] = QP_triangle.get_qp(ii, 2);  // r = t'
        qp[4*ii + 1] = QP_triangle.get_qp(ii, 0);  // s = r'
        qp[4*ii + 2] = QP_triangle.get_qp(ii, 1);  // t = s'
      }
      break;
    
    case 1: // r = 0 : node0 = node0', node3 = node1', node2 = node2'
      for(unsigned int ii{0}; ii < num_pts; ++ii)
      {
        qp[4*ii + 1] = QP_triangle.get_qp(ii, 1);  // s = s'
        qp[4*ii + 2] = QP_triangle.get_qp(ii, 0);  // t = r'
        qp[4*ii + 3] = QP_triangle.get_qp(ii, 2);  // u = t'
      }
      break;
    case 2: // s = 0 : node0 = node0', node1 = node1', node3 = node2'
      for(unsigned int ii{0}; ii < num_pts; ++ii)
      {
        qp[4*ii + 0] = QP_triangle.get_qp(ii, 0);  // r = r'
        qp[4*ii + 2] = QP_triangle.get_qp(ii, 1);  // t = s'
        qp[4*ii + 3] = QP_triangle.get_qp(ii, 2);  // u = t'
      }
      break;
    case 3: // t = 0 : node0 = node0', node2 = node1', node1 = node2'
      for(unsigned int ii{0}; ii < num_pts; ++ii)
      {
        qp[4*ii + 0] = QP_triangle.get_qp(ii, 1);  // r = s'
        qp[4*ii + 1] = QP_triangle.get_qp(ii, 0);  // s = r'
        qp[4*ii + 3] = QP_triangle.get_qp(ii, 2);  // u = t'
      }
      break;
    default:
      SYS_T::print_fatal_if( face_id < 0 || face_id > 3,
        "Error: QuadPts_Gauss_Tet_face, wrong face id input.\n");
      break;
  }

  qw.resize( num_pts );
  for(unsigned int ii{0}; ii < num_pts; ++ii)
    qw[ii] = QP_triangle.get_qw(ii);
}

QuadPts_Gauss_Tet_face::~QuadPts_Gauss_Tet_face()
{
  VEC_T::clean(qp); VEC_T::clean(qw);
}

void QuadPts_Gauss_Tet_face::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== Gauss Points for Tetrahedron on face "<<face_id<< "======="<<std::endl;
  std::cout<<"Num of pt = "<<num_pts<<std::endl;
  std::cout<<"qp.size() = "<<qp.size()<<std::endl;
  std::cout<<"qw.size() = "<<qw.size()<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<std::setprecision(16)<<qp[4*ii]
      <<'\t'<<qp[4*ii+1]<<'\t'<<qp[4*ii+2]
      <<'\t'<<qp[4*ii+3]<<'\t'<<qw[ii]<<'\n';
  std::cout<<"==========================================="<<std::endl;
}