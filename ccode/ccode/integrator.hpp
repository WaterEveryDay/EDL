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
#include <Eigen/Dense>

using Eigen::VectorXd;


std::vector<VectorXd> combine(std::vector<double> v1, std::vector<VectorXd> v2);

//std::vector<std::vector<double>> rk4(std::vector<double> (*func)(double, std::vector<double>), double t0, double step, std::vector<double> x0, bool (*stopcond)(std::vector<double>));

std::vector<VectorXd> rk4(std::function<VectorXd (double, VectorXd)> func, double t0, double step, VectorXd x0, bool (*stopcond)(VectorXd));

#endif /* integrator_hpp */
