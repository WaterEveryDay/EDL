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
#include "atmospheric_model.hpp"
#include "vec_functions.hpp"
#include <math.h>
#include <iostream>

class EOM {
public:
    virtual std::vector<double> dxdt(double t, std::vector<double> state)=0;
};



class PlanarEOM: public EOM {
    // to do: g as a function of altitude
    // to do: better atm model
    // to do: implement L/D
    
    double g0; // gravity constant of the body
    double r0; // radius of the body
    ATMModel *atmmodel; // atmospheric model of the body
    double beta; // ballistic coeff of the vehicle
    double mu_titan = 6.67430e-11*1.3452e23; // TO CHANGE

public:
    PlanarEOM(double g_in, double r0_in, ATMModel* atmmodel_in, double beta_in) {
        g0 = g_in;
        r0 = r0_in;
        atmmodel = atmmodel_in;
        beta = beta_in;
    }
    
    std::vector<double> dxdt(double t, std::vector<double> state) {
        double v = state[0];
        double gamma = state[1];
        double h = state[2];
        
        double g = mu_titan/(r0+h)/(r0+h);
        double rho =  atmmodel->getDensity(h);
        double dvdt = -rho * v * v / 2. / beta - g * sin(gamma);
        double dgammadt = 1./v * (v * v * cos(gamma) / (r0 + h) - g * cos(gamma));
        double dhdt = v * sin(gamma);
        
        return std::vector<double> {dvdt, dgammadt, dhdt};
    }
};


#endif /* EOM_hpp */
