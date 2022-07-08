//
//  atmospheric_model.hpp
//  project
//
//  Created by Sebastien Henry on 7/5/22.
//

#ifndef atmospheric_model_hpp
#define atmospheric_model_hpp
#include <iostream>
#include <stdio.h>

class ATMModel {       // The class
  public:             // Access specifier
    virtual ~ATMModel() {}
    virtual double getDensity(double h)=0;
};

class ExpATMModel: public ATMModel {
    double rho0;
    double H;
public:
    ExpATMModel(double rho0_in, double H_in) {
        rho0 = rho0_in;
        H = H_in;
    }

    double getDensity(double h) {
        double rho = rho0 * exp(-h/H);
        return rho;
    }
};

#endif /* atmospheric_model_hpp */
