//
//  EOM.hpp
//  project
//
//  Created by Sebastien Henry on 7/5/22.
//

#ifndef EOM_hpp
#define EOM_hpp

#include <stdio.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <list>

#include "atmospheric_model.hpp"
#include "vec_functions.hpp"
#include "vehicle.hpp"
#include "Eigen/Dense"
#include "integrator.hpp"

using Eigen::VectorXd;

class EOM {
public:
    virtual VectorXd dxdt(double t, VectorXd state)=0;
    virtual ATMModel* getAtmmodel()=0;
    virtual Vehicle* getVehicle()=0;
    virtual std::vector<VectorXd> propagate(double t0, double step, VectorXd x0, bool (*stopcond)(VectorXd))=0;
};


class PlanarEOM: public EOM {
    // to do: implement L/D
    
    double r0; // radius of the body
    ATMModel *atmmodel; // atmospheric model of the body
    Vehicle *vehicle; // vehicle


public:
    PlanarEOM(ATMModel* atmmodel_in, Vehicle* vehicle_in) {
        atmmodel = atmmodel_in;
        vehicle = vehicle_in;
        r0 = atmmodel->getBodyRadius();
    }
    
    VectorXd dxdt(double t, VectorXd state) {
        double v = state(0);
        double gamma = state(1);
        double h = state(2);
        
        double g = atmmodel->getGravityAcceleration(h);
        double rho =  atmmodel->getDensity(h);
        double lam = atmmodel->getMeanFreePath(h);
        
        double l = vehicle->getRefLength();
        double Kn = lam/l;
        double beta = vehicle->getBeta(Kn);
        
        double dvdt = -rho * v * v / 2. / beta - g * sin(gamma);
        double dgammadt = 1./v * (v * v * cos(gamma) / (r0 + h) - g * cos(gamma));
        double dhdt = v * sin(gamma);
        VectorXd deriv(3);
        deriv << dvdt, dgammadt, dhdt;
        return deriv;
    }
    
    ATMModel* getAtmmodel() {
        return atmmodel;
    }
    
    Vehicle* getVehicle() {
        return vehicle;
    }
    std::vector<VectorXd> propagate(double t0, double step, VectorXd x0, bool (*stopcond)(VectorXd)) {
        auto fp = std::bind(&EOM::dxdt, this, std::placeholders::_1, std::placeholders::_2);
        return rk4(fp, t0, step, x0, stopcond);
    }
};


#endif /* EOM_hpp */
