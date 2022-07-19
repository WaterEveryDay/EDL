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
    virtual double getPressure(double h)=0;
    virtual double getMeanFreePath(double h)=0;
    virtual double getGravityAcceleration(double h)=0;
    virtual double getBodyRadius()=0;
};

class YelleATMModel: public ATMModel {
    // gas constants
    double kb = 1.380649e-23;
    double NA = 6.02214076e23;
    double R  = kb*NA;
    
    // specific to atmosphere on Titan
    double M = 27.81e-3; // unit here is kg/mol // choose between 2.88 and 2.81
    double R_spec = R/M;
    double d_N2 = 364e-12; // kinetic diameter of N2
    
    // titan planet model
    double G = 6.67430e-11; // N m^2/kg^2
    double M_titan = 1.3452e23; // kg
    double r0 = 2574.73e3; // radius Titan (m)
    double mu = G*M_titan; // gravitational param Titan
    
    // data from Yelle
    const static int n_elem = 10;
    double altitudes [n_elem] = {0, 40e3, 80e3, 120e3, 220e3, 300e3, 550e3, 700e3, 900e3, 1300e3};
    double densities [n_elem] = {5.4492, 0.6742, 0.0480, 0.0119, 9.1436e-04, 1.6348e-04, 1.1083e-06, 4.3363e-08, 1.5563e-09, 1.0575e-11};
    double temps[n_elem] = {93.28, 70.66, 123.61, 151.00, 173.21, 178.19, 135.49, 151.42, 175, 175};
    double Fcs[n_elem] = {1.0347, 1.00666, 1.00018, 1.00003, 1, 1, 1, 1, 1, 1};
public:
    YelleATMModel() {}

    double getDensity(double h) {
        double rho = 0;
        double H;
        for (int i = 0; i<n_elem; i++){
            if (h < altitudes[i+1]) {
                H = -(altitudes[i+1]-altitudes[i])/log(densities[i+1]/densities[i]);
                rho = densities[i] * exp(-(h-altitudes[i])/H);
                return rho;
            }
        }
        if (h==altitudes[n_elem-1]) {
            return densities[n_elem-1];
        }
        return rho;
    }
    
    double getTemperature(double h) {
        double T = 175.00;
        double Lhi;
        for (int i = 0; i<n_elem; i++){
            if (h < altitudes[i+1]) {
                Lhi = (temps[i+1]-temps[i])/(altitudes[i+1]-altitudes[i]);
                T = temps[i] + Lhi * (h - altitudes[i]);
                return T;
            }
        }
        return T;
    }
    
    double getPressure(double h) {
       return getTemperature(h)*getDensity(h)*R_spec/getFc(h);
    }

    double getFc(double h) {
        double Fc = 1;
        double slope;
        for (int i = 0; i<n_elem; i++) {
            if (h < altitudes[i+1]) {
                slope = (Fcs[i+1]-Fcs[i])/(altitudes[i+1]-altitudes[i]);
                Fc = Fcs[i] + slope * (h - altitudes[i]);
                return Fc;
            }
        }
        return Fc;
    }
    

    
    double getMeanFreePath(double h) {
        double lambda = kb * getTemperature(h)/(sqrt(2)*M_PI*pow(d_N2,2)*getPressure(h));
        return lambda;
    }
    
    double getGravityAcceleration(double h) {
        return mu/pow((r0+h), 2);
    }
    
    double getBodyRadius() {
        return r0;
    }
    
};

class TitanExpATMModel: public ATMModel {
    // gas constants
    double kb = 1.380649e-23;
    double NA = 6.02214076e23;
    double R  = kb*NA;
    
    // specific to atmosphere on Titan
    double M = 27.81e-3; // unit here is kg/mol // choose between 2.88 and 2.81
    double R_spec = R/M;
    double d_N2 = 364e-12; // kinetic diameter of N2
    
    // titan planet model
    double G = 6.67430e-11; // N m^2/kg^2
    double M_titan = 1.3452e23; // kg
    double r0 = 2574.73e3; // radius Titan (m)
    double mu = G*M_titan; // gravitational param Titan
    
    // data from Yelle
    const static int n_elem = 2;
    double altitudes [n_elem] = {0, 1300e3};
    double densities [n_elem] = {5.4492, 1.0575e-11};
    double temps[n_elem] = {93.28, 175};
    double Fcs[n_elem] = {1.0347, 1};
public:
    TitanExpATMModel() {}

    double getDensity(double h) {
        double rho = 0;
        double H;
        for (int i = 0; i<n_elem; i++){
            if (h < altitudes[i+1]) {
                H = -(altitudes[i+1]-altitudes[i])/log(densities[i+1]/densities[i]);
                rho = densities[i] * exp(-(h-altitudes[i])/H);
                return rho;
            }
        }
        if (h==altitudes[n_elem-1]) {
            return densities[n_elem-1];
        }
        return rho;
    }
    
    double getTemperature(double h) {
        double T = 175.00;
        double Lhi;
        for (int i = 0; i<n_elem; i++){
            if (h < altitudes[i+1]) {
                Lhi = (temps[i+1]-temps[i])/(altitudes[i+1]-altitudes[i]);
                T = temps[i] + Lhi * (h - altitudes[i]);
                return T;
            }
        }
        return T;
    }
    
    double getPressure(double h) {
       return getTemperature(h)*getDensity(h)*R_spec/getFc(h);
    }

    double getFc(double h) {
        double Fc = 1;
        double slope;
        for (int i = 0; i<n_elem; i++) {
            if (h < altitudes[i+1]) {
                slope = (Fcs[i+1]-Fcs[i])/(altitudes[i+1]-altitudes[i]);
                Fc = Fcs[i] + slope * (h - altitudes[i]);
                return Fc;
            }
        }
        return Fc;
    }
    

    
    double getMeanFreePath(double h) {
        double lambda = kb * getTemperature(h)/(sqrt(2)*M_PI*pow(d_N2,2)*getPressure(h));
        return lambda;
    }
    
    double getGravityAcceleration(double h) {
        return mu/pow((r0+h), 2);
    }
    
    double getBodyRadius() {
        return r0;
    }
    
};

//class ExpATMModel: public ATMModel {
//    double rho0;
//    double H;
//    double T;
//public:
//    ExpATMModel(double rho0_in, double H_in, double T_in) {
//        rho0 = rho0_in;
//        H = H_in;
//        T = T_in;
//    }
//
//    double getDensity(double h) {
//        double rho = rho0 * exp(-h/H);
//        return rho;
//    }
//    double getTemperature(double h) {
//        return T;
//    }
//};



#endif /* atmospheric_model_hpp */
