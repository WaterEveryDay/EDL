//
//  functions.cpp
//  project
//
//  Created by Sebastien Henry on 7/5/22.
//

#include <vector>
#include <iostream>


std::vector <double> operator* (double scalar, std::vector <double> v) {
    std::vector<double> v_out(v.size(), 0.);
    for (int i = 0; i < v.size(); i++) {
        v_out[i] = v[i] * scalar;
    }
    return v_out;
}

std::vector <double> operator* (std::vector <double> v, double scalar) {
    return scalar*v;
}

std::vector <double> operator+ (std::vector <double> v1, std::vector <double> v2) {
    std::vector<double> v_out(v1.size(), 0.);
    for (int i = 0; i < v1.size(); i++) {
        v_out[i] = v1[i] + v2[i];
    }
    return v_out;
}

void printvec(std::vector<double> v) {
    for (double i : v) {
        std::cout << i << ", ";
    }
    std::cout << "\n";
}

void printvec(std::vector<std::vector<double>> v) {
    for (std::vector<double> i : v) {
        printvec(i);
    }
}

