//
//  integrator.cpp
//  project
//
//  Created by Sebastien Henry on 7/5/22.
//


#include <vector>
#include <list>
#include "integrator.hpp"

std::vector<VectorXd> combine(std::vector<double> v1, std::vector<VectorXd> v2) {
    std::vector<VectorXd> v3(v1.size());
    for (int i=0; i<v1.size(); i++) {
        VectorXd v1i(1);
        v1i << v1[i];
        VectorXd v2i = v2[i];
        VectorXd v3i(v1i.size() + v2i.size());
        v3i << v1i, v2i;
        v3[i] = v3i;
    }
    return v3;
}

/**
 * Integrator, 4th order Runge-Kutta
 * dynamically reallocating
 * - in:
 * func: the function to integrate f(t, state)
 * t0: time zero
 * step: time step
 * x0: the initial state
 * stopcond: function that sends a boolean of stopping condition (for example, altitude becomes negative)
 */
// std::vector<std::vector<double>> rk4(std::vector<double> (*func)(double, std::vector<double>), double t0, double step, std::vector<double> x0, bool (*stopcond)(std::vector<double>))

std::vector<VectorXd> rk4(std::function<VectorXd (double, VectorXd)> func, double t0, double step, VectorXd x0, bool (*stopcond)(VectorXd)) {
    
    std::vector<VectorXd> x = {};
    std::vector<double> t = {};
    x.push_back(x0);
    t.push_back(t0);
    for (int i = 0; !stopcond(x[i]); i++) {
        VectorXd k1 = func(t[i], x[i]);
        VectorXd k2 = func(t[i]+step/2., x[i] +  step/2. * k1);
        VectorXd k3 = func(t[i]+step/2., x[i] + step/2. * k2);
        VectorXd k4 = func(t[i]+step, x[i] + step * k3);
        
        x.push_back(x[i] + 1./6. * step * (k1 + 2*k2 + 2*k3 + k4));
        t.push_back(t[i] + step);
    }
    return combine(t, x);
}
