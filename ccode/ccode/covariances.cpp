//
//  covariances.cpp
//  ccode
//
//  Created by Sebastien Henry on 7/19/22.
//

#include "covariances.hpp"


MatrixXd propagate_extended(EOM* eom, MatrixXd iniCov_RSW, VectorXd inix_RSW) {
    // inix in RSW
    //
    float h = 1e-6;
    long dim = inix_RSW.size();
    
    MatrixXd A(3,3);
    
    
    for (int i=0; i < dim; i++) {
        VectorXd dx(dim);
        dx.setZero();
        dx[i] = h;
        VectorXd inix_bis = inix_RSW + dx;
        
        
        VectorXd inir_RSW = inix_RSW(seq(0,2));
        VectorXd iniv_RSW = inix_RSW(seq(3,6));
        
        
        VectorXd R(1,0,0);
        VectorXd S(0,1,0);
        
        double c = iniv_RSW.dot(S);
        double s = iniv_RSW.dot(R);
        
        double v = iniv_RSW.norm();
        double gamma = atan2(s, c);
        double h = inir_RSW[0]-eom->getAtmmodel()->getBodyRadius();
        
        
        VectorXd state0_vgh(v, gamma, h);
        EOM->
    }
    
    MatrixXd finalcov;
    return finalcov;
}
