//
//  vehicle.hpp
//  ccode
//
//  Created by Sebastien Henry on 7/8/22.
//

#ifndef vehicle_hpp
#define vehicle_hpp

#include <stdio.h>
#include <math.h>

class Vehicle {
public:
    virtual ~Vehicle() {};
    virtual double getMass()=0;
    virtual double getBeta(double Kn)=0;
    virtual double getRefArea()=0;
    virtual double getRefLength()=0;
};

class Huygens: public Vehicle {
    double Aref = 5.73; // m^2
    double mass = 320; // kg
    double Dmax = 2.7; // m
    
    double CD_FM = 2.1; // free molecular/rarefied flow
    double CD_cont = 1.5; // continuum
    
public:
    Huygens() {}
    
    double getMass() {
        return mass;
    }
    
    double getRefArea() {
        return Aref;
    }
    
    double getRefLength() {
        return Dmax;
    }
    
    double getBeta(double Kn) {
        return getMass()/getRefArea()/getCD(Kn);
    }
    
    double getCD(double Kn) {
        if (Kn > 10) {
            return CD_FM;
        } else if (Kn < 0.001) {
            return CD_cont;
        } else {
            // transitional -> bridge function
            double phi = M_PI * (3./8. + 1./8. * log10(Kn));
            double CD = CD_cont + (CD_FM - CD_cont) * pow(sin(phi), 2);
            return CD;
        }
    }
    
};
#endif /* vehicle_hpp */
