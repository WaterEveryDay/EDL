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
#include "covariances.hpp"
#include <random>
#include <unsupported/Eigen/MatrixFunctions>
#include <string>

bool negative_altitude(VectorXd x) {
    return x(2) < 0;
}

void write_x_to_csv(std::string filename, std::vector<VectorXd> x, std::string columnnames) {
    std::ofstream myFile(filename);
    myFile << columnnames << "\n";
    long jsize = x[0].size();
    long xsize = x.size();
    for (int i=0; i<xsize; i++) {
        VectorXd xi = x[i];
        for (int j=0; j<jsize; j++) {
            myFile << xi(j);
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

int main(int argc, char* argv[])
{
    // TO CHANGE
    int simu = 1; // if 1: canted else: straight down
    
    double step = 0.1; // s
    
    // Monte-Carlo
    bool do_MC = true;
    double step_MC = 1; // s
    int n_MC = 300;
    
    std::string outputDir;
    std::string outputFileName;
    if (argc>1) {
        outputDir =  argv[1]; // "/.../outputDir/"
        outputFileName = argv[2]; // "example.csv"
    } else {
        outputDir = "/Users/sebastienhenry/Downloads/";
        outputFileName = "entry.csv";
    }
    
    ATMModel* atmmod = new YelleATMModel();
    //ATMModel* atmmod = new TitanExpATMModel();
    
    // END TO CHANGE
    
    Vehicle* huygens = new Huygens();
    
    double t0;
    double v_atm;
    double h0;
    double gamma0;
    double psi0;
    double theta0;
    double phi0;
    MatrixXd P_RSW0_LOST(3,3);
    MatrixXd P_RSW0_DLT(3,3);
    
    // Find covariance
    if (simu==1) { // canted 30 deg
        t0 = -268.48;
        v_atm = 6.02852782e3;
        h0 = 1247.68717e3;
        gamma0 = -65.62*M_PI/180;
        psi0 = 259.895 * M_PI/180;
        theta0 = 185.43*M_PI/180;
        phi0 = -8.61*M_PI/180;
        
        P_RSW0_LOST << 2.5715,   -0.2053,    0.2899, // FROM MATLAB
        -0.2053,    0.3604,   -0.0126,
        0.2899,   -0.0126,    0.0942;
        P_RSW0_LOST = 1e4*P_RSW0_LOST;
        
        P_RSW0_DLT << 3.7196,   -0.4312,   -0.0613,
        -0.4312,    0.6614,    0.0207,
        -0.0613,    0.0207,    0.2509;
        P_RSW0_DLT = 1e4*P_RSW0_DLT;
        
    } else { // straight down
        t0 = -268.48;
        v_atm = 6.02852782e3;
        h0 = 1247.68717e3;
        gamma0 = -65.62*M_PI/180;
        psi0 = 270*M_PI/180; // 259.895 * M_PI/180;
        theta0 = 185.43*M_PI/180;
        phi0 = -8.61*M_PI/180;
        
        P_RSW0_LOST << 3.4348, -0.3448, -0.5469,
                    -0.3448, 0.9283, 0.0549,
                    -0.5469, 0.0549, 0.9808; // FROM MATLAB
        P_RSW0_LOST = 1e4*P_RSW0_LOST;
        P_RSW0_DLT << 3.4461, -0.3484, -0.5617,
                    -0.3484, 0.9371, 0.0519,
                    -0.5617, 0.0519, 1.0142;
        P_RSW0_DLT = 1e4*P_RSW0_DLT;
    }
    

    // outputs atmospheric model
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
    
    // CREATE EOM AND PROPAGATE
    NonPlanarEOM* eom = new NonPlanarEOM(atmmod, huygens);
    VectorXd state0(6);
    state0 << v_atm, gamma0, h0, psi0, theta0, phi0;

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<VectorXd> state_integ = eom->propagate(t0, step, state0, negative_altitude);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "duration of integration = " << duration.count() << " ms\n";
    write_x_to_csv(outputDir+outputFileName, state_integ, "t, v, gamma, h, psi, theta, phi");

//    MatrixXd P_RSWf_LOST = propagate_extended_RSW(eom, P_RSW0_LOST, state0, t0, step, negative_altitude);
//    std::cout<< "P_RSWF LOST" << "\n";
//    std::cout<< P_RSWf_LOST << "\n";
    
    MatrixXd P_RSWf_LOST_lonlat = propagate_extended_RSW2State(eom, P_RSW0_LOST, state0, t0, step, negative_altitude);
    std::cout<< "P_GEO_end LOST (h lon lat)" << "\n";
    std::cout<< P_RSWf_LOST_lonlat << "\n";


//    MatrixXd P_RSWf_DLT = propagate_extended_RSW(eom, P_RSW0_DLT, state0, t0, step, negative_altitude);
//    std::cout<< "P_RSWF DLT" << "\n";
//    std::cout<< P_RSWf_DLT << "\n";

    MatrixXd P_RSWf_DLT_lonlat = propagate_extended_RSW2State(eom, P_RSW0_DLT, state0, t0, step, negative_altitude);
    std::cout<< "P_GEO_end DLT (h lon lat)" << "\n";
    std::cout<< P_RSWf_DLT_lonlat << "\n";
    
    // MONTE-CARLO SIMULATION
    if (do_MC==true) {
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0.0, 1.0);
        
        // MC LOST
        std::ofstream MC_LOST(outputDir+"MC_LOST_"+outputFileName);
        MC_LOST << std::setprecision(9);
        MC_LOST << "h, theta, phi\n";
        
        std::cout << "MC LOST\n";
        std::cout << "||||||||||\n";
        double slice = 0;
        double slice_prev = 0;
        std::cout << "|";
        for (int i =0; i< n_MC; i++) {
            VectorXd dx(3);
            dx << distribution(generator), distribution(generator), distribution(generator);
            
            MatrixXd H_RSW2State = sensRSW2State(state0, atmmod->getBodyRadius());
            
            MatrixXd P_state0(3,3);
            VectorXd dx_RSW = P_RSW0_LOST.sqrt() * dx;
            VectorXd dhdthetadphi = H_RSW2State * dx_RSW;
            VectorXd dstate(6);
            dstate.setZero();
            dstate(2) = dhdthetadphi(0);
            dstate(4) = dhdthetadphi(1);
            dstate(5) = dhdthetadphi(2);
            
            VectorXd state0_modified = state0 + dstate;
            std::vector<VectorXd> state_modified = eom->propagate(t0, step_MC, state0_modified, negative_altitude);
            VectorXd state_modified_end = state_modified[state_modified.size()-1](seq(1,6));
            MC_LOST << state_modified_end(2) << ", " << state_modified_end(4) << ", " << state_modified_end(5) << "\n";
            
            slice = double(i)/n_MC*100;
            if (slice > slice_prev + 10) {
                std::cout << "|";
                slice_prev = slice_prev + 10;
            }
        }
        std::cout << "\n";
        MC_LOST.close();
        
        
        // MC DLT
        std::ofstream MC_DLT(outputDir+"MC_DLT_"+outputFileName);
        MC_DLT << std::setprecision(9);
        MC_DLT << "h, theta, phi\n";
        
        std::cout << "MC DLT\n";
        std::cout << "||||||||||\n";
        slice = 0;
        slice_prev = 0;
        std::cout << "|";
        
        for (int i =0; i< n_MC; i++) {
            VectorXd dx(3);
            dx << distribution(generator), distribution(generator), distribution(generator);
            
            MatrixXd H_RSW2State = sensRSW2State(state0, atmmod->getBodyRadius());
            
            MatrixXd P_state0(3,3);
            VectorXd dx_RSW = P_RSW0_DLT.sqrt() * dx;
            VectorXd dhdthetadphi = H_RSW2State * dx_RSW;
            VectorXd dstate(6);
            dstate.setZero();
            dstate(2) = dhdthetadphi(0);
            dstate(4) = dhdthetadphi(1);
            dstate(5) = dhdthetadphi(2);
            
            VectorXd state0_modified = state0 + dstate;
            std::vector<VectorXd> state_modified = eom->propagate(t0, step_MC, state0_modified, negative_altitude);
            VectorXd state_modified_end = state_modified[state_modified.size()-1](seq(1,6));
            MC_DLT << state_modified_end(2) << ", " << state_modified_end(4) << ", " << state_modified_end(5) << "\n";
            
            slice = double(i)/n_MC*100;
            if (slice > slice_prev + 10) {
                std::cout << "|";
                slice_prev = slice_prev + 10;
            }
        }
        MC_DLT.close();
    }
    
    delete atmmod;
    delete huygens;
    
    return 0;
}



// // UNCOMMENT FOR PLANAR EQ.
//bool negative_altitude(VectorXd x) {
//    return x(2) < 0.;
//}
//
//void write_x_to_csv(std::string filename, std::vector<VectorXd> x, std::string columnnames) {
//    std::ofstream myFile(filename);
//    myFile << columnnames << "\n";
//    long jsize = x[0].size();
//    long xsize = x.size();
//    for (int i=0; i<xsize; i++) {
//        VectorXd xi = x[i];
//        for (int j=0; j<jsize; j++) {
//            myFile << xi(j);
//            if (j<jsize-1) {
//                myFile << ", ";
//            }
//            else {
//                myFile << "\n";
//            }
//        }
//    }
//
//    myFile.close();
//}
//
//// Driver Code
//int main(int argc, char* argv[])
//{
//
//    double v_atm = 6.02852782e3;
//    double h0 = 1247.68717e3;
//    double gamma0 = -65.62*M_PI/180;
//
//    double t0 = -268.48;
//    double step = 0.1;
//
//    std::string outputDir;
//    std::string outputFileName;
//    if (argc>1) {
//        outputDir =  argv[1]; // "/.../outputDir/"
//        outputFileName = argv[2]; // "example.csv"
//    } else {
//        outputDir = "/Users/sebastienhenry/Downloads/";
//        outputFileName = "entry.csv";
//    }
//
//    ATMModel* atmmod = new YelleATMModel();
//    //ATMModel* atmmod = new TitanExpATMModel();
//    Vehicle* huygens = new Huygens();
//
//    std::ofstream myAtmFile(outputDir+"atm_model_"+outputFileName);
//    myAtmFile << "h, "<< "rho, " << "T, " << "P, " << "kn, " << "C_D" << "\n";
//    for (double h=0; h<=1300e3; h+=10e3) {
//        myAtmFile << h << ", ";
//        myAtmFile << atmmod->getDensity(h) << ", ";
//        myAtmFile << atmmod->getTemperature(h) << ", ";
//        myAtmFile << atmmod->getPressure(h) << ", ";
//        double lam = atmmod->getMeanFreePath(h);
//        double L = huygens->getRefLength();
//        double Kn = lam/L;
//        myAtmFile << Kn << ", ";
//        double C_D = huygens->getMass() / huygens->getBeta(Kn) / huygens->getRefArea();
//        myAtmFile << C_D << "\n";
//    }
//    myAtmFile.close();
//
//    PlanarEOM eom = PlanarEOM(atmmod, huygens);
//    VectorXd state0(3);
//    state0 << v_atm, gamma0, h0;
//
//
////
//    auto start = std::chrono::high_resolution_clock::now();
//    std::vector<VectorXd> state_integ = eom.propagate(t0, step, state0, negative_altitude);
//    auto stop = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
////
//    std::cout << "duration of integration = " << duration.count() << " ms\n";
//    write_x_to_csv(outputDir+outputFileName, state_integ, "t, v, gamma, h");
//
//    delete atmmod;
//    delete huygens;
//
//    return 0;
//}

