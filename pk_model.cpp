#include "pk_model.h"
#include <iostream>
using namespace std;

// Summary for central compartment:
// Inputs: time vector t, state vectors y, parameters p.
// Returns: PKSummary with Cmax, Tmax, half_life, and AUC for the central compartment.
PKSummary model_summary(const vector<double>& t,const vector<vector<double>>& y,const PKParams& p)
{
    PKSummary s{};
    size_t N = t.size();
    int c_index = (p.n == 3) ? 1 : 0;

    vector<double> C(N);
    C[0] = y[0][c_index];
    s.Cmax = C[0];
    s.Tmax = t[0];

    for (size_t i = 1; i < N; ++i) {
        C[i] = y[i][c_index];
        if (C[i] > s.Cmax) {
            s.Cmax = C[i];
            s.Tmax = t[i];
        }
    }

    s.AUC = trapezoidal_integral(t, C);

    double half_val = 0.5 * s.Cmax;
    s.half_life = 0.0;
    for (size_t i = 0; i < N; ++i)
        if (t[i] > s.Tmax && C[i] <= half_val) {
            s.half_life = t[i] - s.Tmax;
            break;
        }

    return s;
}

// 2-compartment IV bolus model:
// Inputs: parameters p, t_end, h, output stream out.
// Outputs: t and y filled with simulation results; out receives header and CSV lines.
// Effect: sets p.n = 2, runs simulation, writes A1,A2 and C1=A1.
void two_compartment_model(PKParams& p,double t_end,double h,ostream& out,vector<double>& t,vector<vector<double>>& y)
{
    p.n = 2;

    p.kelim = 0.20;
    p.k12   = 0.10;
    p.k21   = 0.05;

    vector<double> y0(2);
    y0[0] = p.dose;   // A1
    y0[1] = 0.0;      // A2
    vector<double> elim_rate;
    
    simulate_rk4(t_end, h, y0, p, t, y, elim_rate);
    //calculate error
    double total_elimination = trapezoidal_integral(t, elim_rate);
    double residual = 0;
    double total_amount = 0.;
    for (double i : y.back()){
        residual += i;
        }
    if (p.dose_steps != 0){//multi-dose/infusion
        total_amount = p.F*(p.dose*((t_end/h)/p.dose_steps + 1.) + p.R_in * t_end + p.R_gut * t_end);
        }
    else { // single dose
        total_amount = p.F*(p.dose + p.R_in * t_end + p.R_gut * t_end);
    }
    double error = total_amount - residual - total_elimination;
    if (error < 0){
      error *= -1.;
      }
    
    PKSummary s = model_summary(t, y, p);
    
    out << "# 3-compartment (gut + central + peripheral)\n"
        << "# Tmax = " << s.Tmax << "hrs" << "\n"
        << " Max amount in central = " << s.Cmax << "mg" << "\n"
        << " half_life = " << s.half_life << "hrs" << "\n"
        << " AUC = " << s.AUC << "mg*hrs" << "\n"
        << " Total dose = " << total_amount << "mg" << "\n" 
        << " residual = " << residual << "mg" << "\n"
        << " total elimination = " << total_elimination << "mg" << "\n"  
        << " error = " << error << "mg" << "\n"
        << " percent error = " << 100. * error/total_amount << "%\n";
    /*out << "time,Ag,A1,A2,C1\n";
    
    /*out << "time,A1,A2,C1\n";

    for (size_t i = 0; i < t.size(); ++i) {
        double A1 = y[i][0];
        double A2 = y[i][1];
        double C1 = A1;
        out << t[i] << "," << A1 << "," << A2 << "," << C1 << "\n";
    }*/
}


// 3-compartment (gut + central + peripheral) model:
// Inputs: parameters p, t_end, h, output stream out.
// Outputs: t and y filled with simulation results; out receives header and CSV lines.
// Effect: sets p.n = 3,runs simulation with Ag,A1,A2, writes C1=A1.
void three_compartment_model(PKParams& p,double t_end, double h, ostream& out, vector<double>& t, vector<vector<double>>& y)
{
    p.n = 3;

    vector<double> y0(3);
    y0[0] = p.F*p.dose;
    y0[1] = 0.0;
    y0[2] = 0.0;
    vector<double> elim_rate;
    
    simulate_rk4(t_end, h, y0, p, t, y, elim_rate);
    
    //calculate error
    double total_elimination = trapezoidal_integral(t, elim_rate);
    double residual = 0.;
    double total_amount = 0.;
    for (double i : y.back()){
        residual += i;
        }
    if (p.dose_steps != 0){//multi-dose/infusion
        total_amount = p.F*(p.dose*((t_end/h)/p.dose_steps + 1.) + p.R_in * t_end + p.R_gut * t_end);
        }
    else { // single dose
        total_amount = p.F*(p.dose + p.R_in * t_end + p.R_gut * t_end);
    }
    double error = total_amount - residual - total_elimination;
    if (error < 0){
      error *= -1.;
      }
    PKSummary s = model_summary(t, y, p);

    out << "# 3-compartment (gut + central + peripheral)\n"
        << "# Tmax = " << s.Tmax << "hrs" << "\n"
        << " Max amount in central = " << s.Cmax << "mg" << "\n"
        << " half_life = " << s.half_life << "hrs" << "\n"
        << " AUC = " << s.AUC << "mg*hrs" << "\n"
        << " total dose = " << total_amount << "mg" << "\n" 
        << " residual = " << residual << "mg" << "\n"
        << " total elimination = " << total_elimination << "mg" << "\n"  
        << " error = " << error << "mg" << "\n"
        << " percent error = " << 100. * error/total_amount << "%\n";

    /*out << "time,Ag,A1,A2,C1\n";

    for (size_t i = 0; i < t.size(); ++i) {
        double Ag = y[i][0];
        double A1 = y[i][1];
        double A2 = y[i][2];
        double C1 = A1;
        out << t[i] << "," << Ag << "," << A1 << "," << A2 << "," << C1 << "\n";
    } */
}
