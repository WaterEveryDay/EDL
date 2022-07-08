//
//  integrator.cpp
//  project
//
//  Created by Sebastien Henry on 7/5/22.
//

#include <vector>
#include "vec_functions.hpp"

std::vector<std::vector<double>> combine(std::vector<double> v1, std::vector<std::vector<double>> v2) {
    std::vector<std::vector<double>> v3; //(v1.size(), std::vector<double>(v1[0].size() + v2[0].size(), 0));
    for (int i=0; i<v1.size(); i++) {
        v3.push_back(std::vector<double> {v1[i]});
        v3[i].insert(v3[i].end(), v2[i].begin(), v2[i].end());
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

std::vector<std::vector<double>> rk4(std::function<std::vector<double>(double, std::vector<double>)> func, double t0, double step, std::vector<double> x0, bool (*stopcond)(std::vector<double>)) {
    
    std::vector<std::vector<double>> x(1, x0);
    std::vector<double> t(1, t0);
    x[0] = x0;
    for (int i = 0; !stopcond(x[i]); i++) {
        std::vector<double> k1 = func(t[i], x[i]);
        std::vector<double> k2 = func(t[i]+step/2., x[i] + step/2. * k1);
        std::vector<double> k3 = func(t[i]+step/2., x[i] + step/2. * k2);
        std::vector<double> k4 = func(t[i]+step, x[i] + step * k3);
        
        x.push_back(x[i] + 1./6. * step * (k1 + 2*k2 + 2*k3 + k4));
        t.push_back(t[i] + step);
    }
    return combine(t,x);
}
