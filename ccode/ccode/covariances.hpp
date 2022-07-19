//
//  covariances.hpp
//  ccode
//
//  Created by Sebastien Henry on 7/19/22.
//

#ifndef covariances_hpp
#define covariances_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include "EOM.hpp"

#endif /* covariances_hpp */

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::seq;

MatrixXd propagate_extended(std::function<std::vector<double>(double, std::vector<double>)> func, MatrixXd covMat, VectorXd iniState);


