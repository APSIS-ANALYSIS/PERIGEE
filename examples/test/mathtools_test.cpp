#include <iostream>
#include <vector> 
#include <random> 
#include <array>

#include "Math_Tools.hpp"

int main( int argc, char * argv[] )
{
    MATH_T::Matrix_Dense<4> Mat1;   // Constructor

    MATH_T::Matrix_SymPos_Dense<4> Mat_SymPos1;   // Constructor 

    MATH_T::Matrix_Dense<4> Mat2(Mat1);      // Default Copy Constructor, No print

    MATH_T::Matrix_SymPos_Dense<4> Mat_SymPos2(Mat2);  // Copy Constructor 

    MATH_T::Matrix_SymPos_Dense<4> Mat_SymPos3(Mat_SymPos1);  // Copy Constructor 

    std::cout << "\nMat1: " << std::endl;
    Mat1.print_info();   

    std::cout << "\nMat2: " << std::endl;
    Mat2.print_info();    

    std::cout << "\nMat_SymPos1: " << std::endl;
    Mat_SymPos1.print_info();  

    std::cout << "\nMat_SymPos2: " << std::endl;
    Mat_SymPos2.print_info();  

    std::cout << "\nMat_SymPos3: " << std::endl;
    Mat_SymPos3.print_info();   

    // ----------------------------

    Mat1.gen_rand();          

    Mat2.gen_rand(); 

    Mat_SymPos1.gen_rand();

    Mat_SymPos2.gen_rand();

    Mat_SymPos3.gen_rand();

    std::cout << "\nMat1: " << std::endl;
    Mat1.print_info();    

    std::cout << "\nMat2: " << std::endl;
    Mat2.print_info();   

    std::cout << "\nMat_SymPos1: " << std::endl;
    Mat_SymPos1.print_info();   

    std::cout << "\nMat_SymPos2: " << std::endl;
    Mat_SymPos2.print_info();   

    std::cout << "\nMat_SymPos3: " << std::endl;
    Mat_SymPos3.print_info();   

    //----------------------------

    std::array<double, 3*3> AA{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};   

    std::array<double, 3*3> BB{1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 7.0};

    MATH_T::Matrix_Dense<3> Mat3(AA);   // Constructor

    MATH_T::Matrix_SymPos_Dense<3> Mat_SymPos4(BB);   // Constructor 

    MATH_T::Matrix_Dense<3> Mat4;      

    MATH_T::Matrix_SymPos_Dense<3> Mat_SymPos5;  

    std::cout << "\nMat3: " << std::endl;
    Mat3.print_info();    

    std::cout << "\nMat4: " << std::endl;
    Mat4.print_info();    

    std::cout << "\nMat_SymPos4: " << std::endl;
    Mat_SymPos4.print_info();   

    std::cout << "\nMat_SymPos5: " << std::endl;
    Mat_SymPos5.print_info();   

    // --------------------------------------

    Mat4 = Mat3;   // Assignment operator =

    std::cout << "\nMat4: " << std::endl;
    Mat4.print_info();    

    // --------------------------------------

    Mat4 = Mat_SymPos4;    // Assigment operator = 

    std::cout << "\nMat4: " << std::endl;
    Mat4.print_info();    

    //---------------------------------------
   
    Mat_SymPos5 = Mat_SymPos4;   // Assigment operator =

    std::cout << "\nMat_SymPos5: " << std::endl;
    Mat_SymPos5.print_info();   

    //-----------------------------------------


    MATH_T::Matrix_Dense<3> Mat5;      

    MATH_T::Matrix_SymPos_Dense<3> Mat_SymPos6; 

    MATH_T::Matrix_Dense<3> Mat6; 

    MATH_T::Matrix_Dense<3> Mat, MM, NN; 


    for (int ii = 0; ii < 1; ++ii)
    {
        Mat5.gen_rand();
        Mat.gen_rand();
        Mat_SymPos6 = Mat5;   // Copy Constructor -> Assignment operator =

        std::cout << "\nMat_SymPos6: " << std::endl;
        Mat_SymPos6.print_info();   

        std::cout << "\nMat5: " << std::endl;
        Mat5.print_info();   

        std::cout << "\nMat: " << std::endl;
        Mat.print_info(); 

        // -----------------------------------

        Mat6 = Mat_SymPos6;

        MM.Mult(Mat5, Mat);
        NN.Mult(Mat6, Mat);

        std::cout << "\nMat6: " << std::endl;
        Mat6.print_info();   

        std::cout << "\nMM: " << std::endl;
        MM.print_info();  

        std::cout << "\nNN: " << std::endl;
        NN.print_info();  

        // ----------------------------------

        for (int jj = 0; jj < 3 * 3; ++jj)
        {
            if( !MATH_T::equals( MM(jj), NN(jj), 1.0e-15) )
            {
                std::cout<<"error: Matrix MM("<<jj<<") does not match NN("<<jj<<"). \n";
            }
        }
    }


    return 0;


}

// EOF