#ifndef QUADPTS_GAUSS_TET_FACE_HPP
#define QUADPTS_GAUSS_TET_FACE_HPP
// ==================================================================
// QuadPts_Gauss_Tet.hpp
// The Gaussian quadrature rule for a triangular domain expressed on
// a face of tet element defined by four vertex points:
// [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]
//
//                     t
//                     ^
//                     |
//                     3
//                    /| `.
//                   / |    `.
//                  /  |       `.
//                 /   |          `.
//                /    |             `.     
//               /     |                `.
//              /      |                   `.   
//             /      ,0 - - - - - - - - - - -`2 - - -> s
//            /     ,'  (u)               ,  "
//           /    ,'                ,  "
//          /   ,'            ,  "
//         /  ,'        ,  "
//        / ,'    ,  "
//       /,',  "
//      1'
//    ,'
// r * 
//  
// Face 0:  u = 0
//      t                                        s'
//      ^                                        ^
//      |                                        |
//      3                                        2'
//      |  `.                      map           |  `.
//      |     `.                  <----          |     `.  
//      |        `.                              |        `.
//      |           `.                           |           `.
//      |              `.                        |              `.
//  (r) 1 - - - - - - - - 2 - - -> s        (t') 0'- - - - - - - - 1'- - -> r'
//
// Face 1:  r = 0
//      s                                        s'
//      ^                                        ^
//      |                                        |
//      2                                        2'
//      |  `.                      map           |  `.
//      |     `.                  <----          |     `.  
//      |        `.                              |        `.
//      |           `.                           |           `.
//      |              `.                        |              `.
//  (u) 0 - - - - - - - - 3 - - -> t        (t') 0'- - - - - - - - 1'- - -> r'
//
// Face 2:  s = 0
//      t                                        s'
//      ^                                        ^
//      |                                        |
//      3                                        2'
//      |  `.                      map           |  `.
//      |     `.                  <----          |     `.  
//      |        `.                              |        `.
//      |           `.                           |           `.
//      |              `.                        |              `.
//  (u) 0 - - - - - - - - 1 - - -> r        (t') 0'- - - - - - - - 1'- - -> r'
//
// Face 3:  t = 0
//      r                                        s'
//      ^                                        ^
//      |                                        |
//      1                                        2'
//      |  `.                      map           |  `.
//      |     `.                  <----          |     `.  
//      |        `.                              |        `.
//      |           `.                           |           `.
//      |              `.                        |              `.
//  (u) 0 - - - - - - - - 2 - - -> s        (t') 0'- - - - - - - - 1'- - -> r'
//
// Date Created: Oct. 24 2023
// ==================================================================
#include "QuadPts_Gauss_Triangle.hpp"

class QuadPts_Gauss_Tet_face : public IQuadPts
{
    public:
      QuadPts_Gauss_Tet_face( const int &in_num_pts, const int & in_face_id );

      virtual ~QuadPts_Gauss_Tet_face();

      virtual void print_info() const;

      virtual int get_dim() const {return 4;}

      virtual int get_boundary_id() const {return face_id;}

      virtual int get_num_quadPts() const {return num_pts;}

      virtual double get_qp(unsigned int ii, unsigned int comp) const
      {return qp[4*ii+comp];}

      virtual double get_qw(unsigned int ii) const
      {return qw[ii];}

      virtual IQuadPts * get_lower_QP(){return &QP_triangle;}

    private:
        const int num_pts;

        // qp : length 4 * num_pts. Stores the r-s-t-u coordinates of the 
        //      quadrature points.
        //      u = 1 - r - s - t
        // qw : length num_pts. Stores the quadrature weights.
        std::vector<double> qp {};
        std::vector<double> qw {};

        const int face_id;

        QuadPts_Gauss_Triangle QP_triangle;
};

#endif