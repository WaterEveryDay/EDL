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
#include <math.h>

#endif /* covariances_hpp */

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::seq;

MatrixXd propagate_extended_RSW(EOM* eom, MatrixXd cov_init_RSW, VectorXd state_init, double t0, double step, bool (*stopcond)(VectorXd));
MatrixXd propagate_extended_state(EOM* eom, MatrixXd cov_init_state, VectorXd state_init, double t0, double step, bool (*stopcond)(VectorXd));
MatrixXd propagate_extended_RSW2State(EOM* eom, MatrixXd cov_init_RSW, VectorXd state_init, double t0, double step, bool (*stopcond)(VectorXd));
//VectorXd state2Titan(VectorXd state, double r0);
//VectorXd state2RSW(VectorXd state, double r0);
//MatrixXd Titan2RSW(VectorXd state, double r0);

MatrixXd sensState2RSW(VectorXd state, double r0);
MatrixXd sensRSW2State(VectorXd state, double r0);





