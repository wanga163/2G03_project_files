#ifndef PK_MODEL_H
#define PK_MODEL_H

#include <vector>
#include <ostream>
using namespace std;

//structs
struct PKParams {
    int n; // number of compartments (2 or 3)
    double dose; // bolus dose
    int dose_steps; // steps between repeated doses (0 = no repeats)
    int dose_compartment; // index to receive dose
    double F; // bioavailability (for oral)

    // rate constants
    double kabsorp; // gut -> central
    double kelim; // elimination from central
    double k12; // central -> peripheral
    double k21; // peripheral -> central
    double R_in; // IV infusion rate (mg/h)
    double R_gut; // zero-order oral release rate (mg/h)
    
    //Effective dose calculation
    double Vd; // effective volume of distribution (L) (224L for 70 kg)
    double Ceffective; // effective plasma concentration (mg/L) (0.03mg/L)

    PKParams();
};

struct PKSummary {
    double Cmax;
    double Tmax;
    double half_life;
    double AUC;
};

// Calculus (functions.cpp)
double trapezoidal_integral(const vector<double>& x, const vector<double>& y);

void pk_derivatives(double t, const vector<double>& y,vector<double>& dydt, const PKParams& p);

void rk4_step(double t, double h, vector<double>& y, const PKParams& p);

void simulate_rk4(double t_end, double h,const vector<double>& y0,PKParams& p,vector<double>& t,vector<vector<double>>& y, vector<double>& elim_rate);


// 2/3 compartment model + summary function (pk_model.cpp)
PKSummary model_summary(const vector<double>& t, const vector<vector<double>>& y,const PKParams& p);

void two_compartment_model(PKParams& p,double t_end,double h,ostream& out,vector<double>& t,vector<vector<double>>& y);

void three_compartment_model(PKParams& p,double t_end,double h,ostream& out,vector<double>& t,vector<vector<double>>& y);

#endif // PK_MODEL_H
