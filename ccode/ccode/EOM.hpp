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
#include "atmospheric_model.hpp"
#include "vec_functions.hpp"
#include "vehicle.hpp"

class EOM {
public:
    virtual std::vector<double> dxdt(double t, std::vector<double> state)=0;
    ATMModel* getAtmmodel()=0;
    
    Vehicle* getVehicle()=0;
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
    
    std::vector<double> dxdt(double t, std::vector<double> state) {
        double v = state[0];
        double gamma = state[1];
        double h = state[2];
        
        double g = atmmodel->getGravityAcceleration(h);
        double rho =  atmmodel->getDensity(h);
        double lam = atmmodel->getMeanFreePath(h);
        
        double l = vehicle->getRefLength();
        double Kn = lam/l;
        double beta = vehicle->getBeta(Kn);
        
        double dvdt = -rho * v * v / 2. / beta - g * sin(gamma);
        double dgammadt = 1./v * (v * v * cos(gamma) / (r0 + h) - g * cos(gamma));
        double dhdt = v * sin(gamma);
        
        return std::vector<double> {dvdt, dgammadt, dhdt};
    }
    
    ATMModel* getAtmmodel() {
        return atmmodel;
    }
    
    Vehicle* getVehicle() {
        return vehicle;
    }
    
};


#endif /* EOM_hpp */
