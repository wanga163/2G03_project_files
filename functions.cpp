#include "pk_model.h"
#include <vector>
using namespace std;
// PKParams constructor:
// Initializes default values; rate constants set manually.
PKParams::PKParams()
    : n(2),
      dose(10.0),
      dose_steps(0),
      dose_compartment(0),
      F(0.9),
      kabsorp(0.0),
      kelim(0.0),
      k12(0.0),
      k21(0.0),
      R_in(0.0),
      R_gut(0.0),
      Vd(224.0),
      Ceffective(0.03)
{
}


// Trapezoidal rule:
// Inputs: x (time points), y (values at those times).
// Returns: scalar numerical integral of y with respect to x.
double trapezoidal_integral(const vector<double>& x,const vector<double>& y)
{
    size_t n = x.size();

    double A = 0.0;
    for (size_t i = 0; i + 1 < n; ++i)
        A += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]); // (x1+x2) * height/2

    return A;
}

// ODE system:
// Inputs: time t 9not used), current state y, parameters p.
// Output: dydt is resized and overwritten with dy/dt for each component of y.
void pk_derivatives(double,const vector<double>& y,vector<double>& dydt,const PKParams& p)
{
    size_t n = y.size();
    dydt.assign(n, 0.0);

    if (p.n == 2) {
        double A1 = y[0], A2 = y[1];
        double input = p.R_in;
        dydt[0] = input -(p.kelim + p.k12) * A1 + p.k21 * A2;
        dydt[1] =  p.k12 * A1 - p.k21 * A2;
    }
    else if (p.n == 3) {
        double Ag = y[0], A1 = y[1], A2 = y[2];
        double input = p.F*p.R_gut;
        double absorption = p.kabsorp * Ag;
        double elimination = p.kelim * A1;

        dydt[0] = input - absorption;
        dydt[1] = p.F * absorption - elimination - p.k12 * A1 + p.k21 * A2;
        dydt[2] = p.k12 * A1 - p.k21 * A2;
    }
}

// RK4 step:
// Inputs: current time t, step size h, state vector y, parameters p.
// Effect: y is updated in place to y(t + h) using one RK4 step.
void rk4_step(double t, double h,vector<double>& y,const PKParams& p)
{
    size_t n = y.size();
    vector<double> k1(n), k2(n), k3(n), k4(n), temp(n);

    pk_derivatives(t, y, k1, p);

    for (size_t i = 0; i < n; ++i)
        temp[i] = y[i] + 0.5 * h * k1[i]; 
    pk_derivatives(t + 0.5*h, temp, k2, p);

    for (size_t i = 0; i < n; ++i)
        temp[i] = y[i] + 0.5 * h * k2[i];
    pk_derivatives(t + 0.5*h, temp, k3, p);

    for (size_t i = 0; i < n; ++i)
        temp[i] = y[i] + h * k3[i];
    pk_derivatives(t + h, temp, k4, p);

    for (size_t i = 0; i < n; ++i)
        y[i] += (h/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}

// Time integration with RK4:
// Inputs: t_end, step size h, initial state y0, parameters p.
// Outputs: t is filled with time points, y is filled with state vectors at each time.
// y and t are resized and filled; y dosing applied to yc.
void simulate_rk4(double t_end, double h,const vector<double>& y0,PKParams& p, vector<double>& t,vector<vector<double>>& y, vector<double>& elim_rate)
{
    int N = static_cast<int>(t_end / h) + 1;
    size_t n = y0.size();

    t.resize(N);
    y.assign(N, vector<double>(n));
    elim_rate.resize(N);
    
    double time = 0.0;
    vector<double> yc = y0;
    t[0] = 0.0;
    y[0] = yc;
    
    int c_index = 0;
    if (p.n == 3) c_index = 1; //position of central in vector
    elim_rate[0] = p.kelim * yc[c_index];
    
    for (int i = 1; i < N; ++i) {
        rk4_step(time, h, yc, p);
        time += h;

        if (p.dose_steps > 0 &&
            (i % p.dose_steps == 0) && p.dose_compartment >= 0 && p.dose_compartment < static_cast<int>(n)) yc[p.dose_compartment] += p.dose; //add multiple injection

        t[i] = time;
        y[i] = yc;
        elim_rate[i] = p.kelim * yc[c_index];
    }
}
