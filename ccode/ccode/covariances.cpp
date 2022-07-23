//
//  covariances.cpp
//  ccode
//
//  Created by Sebastien Henry on 7/19/22.
//

#include "covariances.hpp"

MatrixXd propagate_extended_state(EOM* eom, MatrixXd cov_init_state, VectorXd state_init, double t0, double step, bool (*stopcond)(VectorXd)) {

    double r0 = eom->getAtmmodel()->getBodyRadius();
    
    MatrixXd F(3, 3);
    F.setZero();
    
    MatrixXd H_RSW2state_init = sensRSW2State(state_init, r0);
    std::vector<VectorXd> states = eom->propagate(t0, step, state_init, stopcond);
    VectorXd state_end = states[states.size()-1](seq(1,6))  ;
    
    float dstep_dist = 1000; // m
    float dstep_angle = 1*M_PI/180; // rad
    
    for (int i=0; i < 3; i++) {

        VectorXd dhdthetadphi(3);
        dhdthetadphi.setZero();
        if (i==0) {
            dhdthetadphi(i) = dstep_dist;
        } else {
            dhdthetadphi(i) = dstep_angle;
        }
        
        VectorXd dstate_init(6);
        dstate_init.setZero();
        dstate_init(2) = dhdthetadphi(0);
        dstate_init(4) = dhdthetadphi(1);
        dstate_init(5) = dhdthetadphi(2);
        
        VectorXd state_init_modified = state_init + dstate_init;
        std::vector<VectorXd> states_modified = eom->propagate(t0, step, state_init_modified, stopcond);
        VectorXd state_end_modified = states_modified[states_modified.size()-1](seq(1,6));
        VectorXd dstate_final = state_end_modified - state_end;
        
        double dhdi = dstate_final(2)/dhdthetadphi(i);
        double dthetadi = dstate_final(4)/dhdthetadphi(i);
        double dphidi = dstate_final(5)/dhdthetadphi(i);
        
        F(0,i) = dhdi;
        F(1,i) = dthetadi;
        F(2,i) = dphidi;
    }

    MatrixXd cov_end_RSW = F * cov_init_state * F.transpose();
    return cov_end_RSW;
}

MatrixXd propagate_extended_RSW2State(EOM* eom, MatrixXd cov_init_RSW, VectorXd state_init, double t0, double step, bool (*stopcond)(VectorXd)) {

    double r0 = eom->getAtmmodel()->getBodyRadius();
    
    MatrixXd H_RSW2state_init = sensRSW2State(state_init, r0);
    MatrixXd cov_init_state =  H_RSW2state_init * cov_init_RSW * H_RSW2state_init.transpose();
    MatrixXd cov_end_state = propagate_extended_state(eom, cov_init_state, state_init, t0, step, stopcond);
    return cov_end_state;
}

MatrixXd propagate_extended_RSW(EOM* eom, MatrixXd cov_init_RSW, VectorXd state_init, double t0, double step, bool (*stopcond)(VectorXd)) {

    double r0 = eom->getAtmmodel()->getBodyRadius();
    
    std::vector<VectorXd> states = eom->propagate(t0, step, state_init, stopcond);
    VectorXd state_end = states[states.size()-1](seq(1,6))  ;
    MatrixXd H_state2RSW_end = sensState2RSW(state_end, r0);
    
    MatrixXd cov_end_state = propagate_extended_RSW2State(eom, cov_init_RSW, state_init, t0, step, stopcond);
    MatrixXd cov_end_RSW = H_state2RSW_end * cov_end_state * H_state2RSW_end.transpose();
    return cov_end_RSW;
}


MatrixXd sensRSW2State(VectorXd state, double r0){
    // double v = state(0); // speed
    // double gamma = state(1); // flight path angle
    double h = state(2); // altitude
    double psi = state(3); // heading angle
    // double theta = state(4); // longitude
    double phi = state(5); // latitude
    double r = h+r0;
    MatrixXd A(3,3);
    double dhdR = 1;
    double dhdS = 0;
    double dhdW = 0;
    
    double dthetadR = 0;
    double dthetadS = sin(-psi)/(r*cos(phi));
    double dthetadW = sin(-psi+M_PI/2) / (r*cos(phi));
    
    double dphidR = 0;
    double dphidS = cos(psi)/r;
    double dphidW = cos(psi-M_PI/2)/r;
    
    A << dhdR, dhdS, dhdW, dthetadR, dthetadS, dthetadW, dphidR, dphidS, dphidW;
    
    return A;
}

MatrixXd sensState2RSW(VectorXd state, double r0) {
    // double v = state(0); // speed
    // double gamma = state(1); // flight path angle
    double h = state(2); // altitude
    double psi = state(3); // heading angle
    // double theta = state(4); // longitude
    double phi = state(5); // latitude
    double r = h+r0;

    MatrixXd A(3,3);

    double dRdh = 1;
    double dRdtheta = 0;
    double dRdphi = 0;

    double dSdh = 0;
    double dSdtheta = r*cos(phi)/sin(-psi);
    double dSdphi = r/cos(psi);

    double dWdh = 0;
    double dWdtheta = r*cos(phi)/sin(-psi+M_PI/2);
    double dWdphi = r/cos(psi-M_PI/2);

    A << dRdh, dRdtheta, dRdphi, dSdh, dSdtheta, dSdphi, dWdh, dWdtheta, dWdphi;
    // MatrixXd A = sensRSW2State(state, r0);
    return A;
}
