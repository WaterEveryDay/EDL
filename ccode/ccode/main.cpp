//  main.cpp
//
//  Created by Sebastien Henry on 7/5/22.
//

#include <iostream>
#include <vector>
#include <math.h>
#include <functional>
#include <chrono>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>

#include "integrator.hpp"
#include "EOM.hpp"
#include "atmospheric_model.hpp"
#include "vehicle.hpp"


bool negative_altitude(std::vector<double> x) {
    return x[2] < 0.;
}

void write_x_to_csv(std::string filename, std::vector<std::vector<double>> x, std::string columnnames) {
    std::ofstream myFile(filename);
    myFile << columnnames << "\n";
    long jsize = x[0].size();
    long xsize = x.size();
    for (int i=0; i<xsize; i++) {
        std::vector<double> xi = x[i];
        for (int j=0; j<jsize; j++) {
            myFile << xi[j];
            if (j<jsize-1) {
                myFile << ", ";
            }
            else {
                myFile << "\n";
            }
        }
    }
    
    myFile.close();
}

// Driver Code
int main(int argc, char* argv[])
{
    
    double v_atm = 6.02852782e3;
    double h0 = 1247.68717e3;
    double gamma0 = -65.62*M_PI/180;
    
    double t0 = -268.48;
    double step = 0.1;
    
    std::string outputDir;
    std::string outputFileName;
    if (argc>1) {
        outputDir =  argv[1]; // "/.../outputDir/"
        outputFileName = argv[2]; // "example.csv"
    } else {
        outputDir = "/Users/sebastienhenry/Downloads/";
        outputFileName = "entry.csv";
    }
    
    //ATMModel* atmmod = new YelleATMModel();
    ATMModel* atmmod = new TitanExpATMModel();
    Vehicle* huygens = new Huygens();
    
    std::ofstream myAtmFile(outputDir+"atm_model_"+outputFileName);
    myAtmFile << "h, "<< "rho, " << "T, " << "P, " << "kn, " << "C_D" << "\n";
    for (double h=0; h<=1300e3; h+=10e3) {
        myAtmFile << h << ", ";
        myAtmFile << atmmod->getDensity(h) << ", ";
        myAtmFile << atmmod->getTemperature(h) << ", ";
        myAtmFile << atmmod->getPressure(h) << ", ";
        double lam = atmmod->getMeanFreePath(h);
        double L = huygens->getRefLength();
        double Kn = lam/L;
        myAtmFile << Kn << ", ";
        double C_D = huygens->getMass() / huygens->getBeta(Kn) / huygens->getRefArea();
        myAtmFile << C_D << "\n";
    }
    myAtmFile.close();
    
    PlanarEOM eom = PlanarEOM(atmmod, huygens);
    std::vector<double> state0 {v_atm, gamma0, h0};

//    auto fp = std::bind(&EOM::dxdt, eom, std::placeholders::_1, std::placeholders::_2);
//
//    auto start = std::chrono::high_resolution_clock::now();
//    std::vector<std::vector<double>> state_integ = rk4(fp, t0, step, {v_atm, gamma0, h0}, negative_altitude);
//    auto stop = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
//
//    std::cout << "duration of integration = " << duration.count() << " ms\n";
    //write_x_to_csv(outputDir+outputFileName, state_integ, "t, v, gamma, h");
    
    delete atmmod;
    delete huygens;
    
    return 0;
}
