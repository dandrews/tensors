//
//  multi_array_demo.cpp
//  
//
//  Created by James Li on 7/19/13.
//
//

#include "multi_array_demo.h"
#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <iostream>
#include <cassert>

int main (){
    // Create a 3D array that is 3x4x2x6
    typedef boost::multi_array<double, 3>array_type;
    typedef array_type::index index;
    boost::array<array_type::index, 3> shape = {{ 3, 4, 2 }};
    array_type A(shape);
    
//    // Assign values to the elements
//    int values =0;
//    for (index i = 0; i != 3; ++i){
//        for (index j = 0; j !=4; ++j){
//            for (index k = 0; k != 2; ++k){
//                for(index l = 0; l != 6; ++l){
//                A[i][j][k][l] = values++;
//                }
//            }  
//        }
//    }
//
//    //     Verify and print values to console
//    int verify = 0;
//    for (index i = 0; i != 3; ++i){
//        for (index j = 0; j!= 4; ++j){
//            for (index  k = 0; k != 2; ++k){
//                for (index l = 0; l != 6; ++l){
//                    assert(A[i][j][k][l]==verify++);
//                   // std::cout << A[i][j][k][l] << " ";
//                }
//            }
//        }
//    }
    
    return 0;
}
