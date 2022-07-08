//
//  integrator.hpp
//  project
//
//  Created by Sebastien Henry on 7/5/22.
//

#ifndef integrator_hpp
#define integrator_hpp

#include <stdio.h>
#include <vector>

std::vector<std::vector<double>> combine(std::vector<double> v1, std::vector<std::vector<double>> v2);

//std::vector<std::vector<double>> rk4(std::vector<double> (*func)(double, std::vector<double>), double t0, double step, std::vector<double> x0, bool (*stopcond)(std::vector<double>));

std::vector<std::vector<double>> rk4(std::function<std::vector<double>(double, std::vector<double>)> func, double t0, double step, std::vector<double> x0, bool (*stopcond)(std::vector<double>));

#endif /* integrator_hpp */
