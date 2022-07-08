//
//  functions.hpp
//  project
//
//  Created by Sebastien Henry on 7/5/22.
//

#ifndef vec_functions_hpp
#define vec_functions_hpp

#include <stdio.h>
#include <vector>

std::vector <double> operator* (double scalar, std::vector <double> v);
std::vector <double> operator* (std::vector <double> v, double scalar);
std::vector <double> operator+ (std::vector <double> v1, std::vector <double> v2);
void printvec(std::vector<double> v);
void printvec(std::vector<std::vector<double>> v);


#endif /* functions_hpp */
