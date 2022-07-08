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

#include "integrator.hpp"
#include "EOM.hpp"
#include "atmospheric_model.hpp"
//#include "print_hi.hpp"

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
    std::string outputDir;
    std::string outputFileName;
    if (argc>1) {
        outputDir =  argv[1]; // "/.../outputDir/"
        outputFileName = argv[2]; // "example.csv"
    } else {
        outputDir = "/Users/sebastienhenry/Downloads/";
        outputFileName = "entry.csv";
    }
    // characteristics of Titan
    double g0 = 1.352; // m/s^2
    double r0 = 2574.73e3; // m
    

    double beta = 28.75; // kg/m^2
    double v_atm = 6.0312e3;
    double h0 = 1270.01e3;
    double gamma0 = -65*M_PI/180;
    
    double t0 = 0;
    double step = 0.1;
    
    // very course model of Titan Atm
    // pressure is around 1.5 Earth one, and T around /3 Earth one
    // rho_titan = 1.5 * 3 * rho_earth around 5.5
    // rho at 1400 is around 1e-12
    // H = -1400e3/log(1e-12/5.5)
    double rho0 = 5.5125;
    double H = 47.7196e3;
    
    ATMModel* titanexp = new ExpATMModel(rho0, H);
    PlanarEOM eom = PlanarEOM(g0, r0, titanexp, beta);
    std::vector<double> state0 {v_atm, gamma0, h0};
    
    auto fp = std::bind(&EOM::dxdt, eom, std::placeholders::_1, std::placeholders::_2);

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double>> state_integ = rk4(fp, t0, step, {v_atm, gamma0, h0}, negative_altitude);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "duration of integration = " << duration.count() << " ms\n";

    write_x_to_csv(outputDir+outputFileName, state_integ, "t, v, gamma, h");
    
    delete titanexp;
    
    return 0;
}
