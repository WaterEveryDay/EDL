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
    virtual double getTemperature(double h)=0;
};


class ExpATMModel: public ATMModel {
    double rho0;
    double H;
    double T;
public:
    ExpATMModel(double rho0_in, double H_in, double T_in) {
        rho0 = rho0_in;
        H = H_in;
        T = T_in;
    }

    double getDensity(double h) {
        double rho = rho0 * exp(-h/H);
        return rho;
    }
    double getTemperature(double h) {
        return T;
    }
};

class YelleATMModel: public ATMModel {
    double altitudes [10] = {0, 40e3, 80e3, 120e3, 220e3, 300e3, 550e3, 700e3, 900e3, 1300e3};
    double densities [10] = {5.4629, 0.6759, 0.0481, 0.0119, 9.1666e-04, 1.6389e-04, 1.1111e-06, 4.3472e-08, 1.5602e-09, 1.0602e-11};
    double temps[10] = {93.28, 70.66, 123.61, 151.00, 173.21, 178.19, 135.49, 151.42, 175, 175};
    
public:
    YelleATMModel() { }

    double getDensity(double h) {
        double rho = 0;
        double H;
        for (int i = 0; i<sizeof(altitudes); i++){
            if (h < altitudes[i+1]) {
                H = -(altitudes[i+1]-altitudes[i])/log(densities[i+1]/densities[i]);
                rho = densities[i] * exp(-(h-altitudes[i])/H);
                return rho;
            }
        }
        return rho;
    }
    
    double getTemperature(double h) {
        double T = 175.00;
        double Lhi;
        for (int i = 0; i<sizeof(altitudes); i++){
            if (h < altitudes[i+1]) {
                Lhi = (temps[i+1]-temps[i])/(altitudes[i+1]-altitudes[i]);
                T = temps[i] + Lhi * (h - altitudes[i]);
                return T;
            }
        }
        return T;
    }
};


#endif /* atmospheric_model_hpp */
