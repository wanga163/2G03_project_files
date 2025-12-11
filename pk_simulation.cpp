#include "pk_model.h"
#include <vector>
#include <iostream>
#include "cpgplot.h"
using namespace std;

// Plot a PK
// t  = time points
// y  = state vectors at each time
// p  = parameters
void plot_run(const vector<double>& t, const vector<vector<double>>& y, const PKParams& p)
{
    int n_points = static_cast<int>(t.size()); //points

    //select compartments
    vector<int> idx;
    vector<const char*> labels;

    if (p.n == 2) {
        idx.push_back(0); labels.push_back("Central");
        idx.push_back(1); labels.push_back("Peripheral");
    } else if (p.n == 3) {
        idx.push_back(0); labels.push_back("Gut");
        idx.push_back(1); labels.push_back("Central");
        idx.push_back(2); labels.push_back("Peripheral");
    } else {
        cerr << "Unexpected number of compartments: " << p.n << "\n";
        return;
    }

    // float conversion
    vector<float> xt(n_points);
    for (int i = 0; i < n_points; ++i) {
        xt[i] = static_cast<float>(t[i]);
    }

    float ymin = 0.0f;
    float ymax = 0.0f;
    bool first = true;

    for (int j : idx) {
      for (int i = 0; i < n_points; ++i) {
          float val = static_cast<float>(y[i][j]);
          if (first) {
              ymin = ymax = val;
              first = false;
          } else {
              if (val < ymin) ymin = val;
              if (val > ymax) ymax = val;
          }
      }
    }
    
    /*double A_effective = p.Ceffective * p.Vd;
    float  eff_line    = static_cast<float>(A_effective);

    if (eff_line < ymin) ymin = eff_line;
    if (eff_line > ymax) ymax = eff_line;*/
    
    float xmin = xt.front();
    float xmax = xt.back();
    float dy   = (ymax - ymin) * 0.1f;
    if (dy <= 0.0f) dy = 1.0f;
    ymin -= dy;
    ymax += dy;

    if (cpgopen(const_cast<char*>("/XWINDOW")) <= 0) {
      cerr << "Unable to open PGPLOT device.\n";
      return;
    }

    cpgenv(xmin, xmax, ymin, ymax, 0, 0);
    cpglab("Time (h)", "Amount (mg)", " ");

    for (size_t k = 0; k < idx.size(); ++k) {
      int j = idx[k];

      vector<float> yt(n_points);
      for (int i = 0; i < n_points; ++i) {
          yt[i] = static_cast<float>(y[i][j]);
      }

      // Color index: 1 = white, 2 = red, 3 = green, etc.
      cpgsci(static_cast<int>(k) + 2);
      cpgline(n_points, xt.data(), yt.data());
    }
    /*
    cpgsci(1); // white
    float hx[2] = { xmin, xmax };
    float hy[2] = { eff_line, eff_line };
    cpgline(2, hx, hy); */

    cpgsci(1);
    float x_legend = xmin + 0.05f * (xmax - xmin);
    float y_legend = ymax - 0.05f * (ymax - ymin);

    for (size_t k = 0; k < labels.size(); ++k) {
      cpgsci(static_cast<int>(k) + 2);
      cpgptxt(x_legend,y_legend - static_cast<float>(k) * 0.05f * (ymax - ymin),0.0f, 0.0f,labels[k]);
    }
    
    /*cpgsci(1);
    cpgptxt(x_legend,y_legend - static_cast<float>(labels.size()) * 0.05f * (ymax - ymin),0.0f, 0.0f,"Effective Conc. Line");
    cpgsci(1);*/
    cpgclos();
}


// IV bolus, 2-compartment model
int run_iv_bolus()
{
    PKParams p;

    p.dose = 10.0; // mg
    p.dose_steps = 0;  // single bolus
    p.dose_compartment = 0; // central for 2-comp
    p.F = 1.0; // not used for IV

    p.kabsorp = 1.16; // no gut compartment in 2-comp IV <---------------------------------------------------------------------------------------- Change these
    p.kelim = 0.44; // elimination from central
    p.k12 = 0.3; // central -> peripheral
    p.k21 = 0.1; // peripheral -> central

    double t_end = 24.0; // hours
    double h = 0.001; // step size

    vector<double> t;
    vector<vector<double>> y;

    two_compartment_model(p, t_end, h, cout, t, y);
    plot_run(t, y, p);

    return 0;
}

// Oral bolus, 3-compartment model (gut + central + peripheral)
int run_oral_bolus()
{
    PKParams p;
    p.dose = 10.0;   // mg
    p.dose_steps = 0;      // single dose 
    p.dose_compartment = 0;      // gut for 3-comp
    p.F = 0.239;    // oral bioavailability

    p.kabsorp = 1.16;   // gut -> central <---------------------------------------------------------------------------------------- Change these
    p.kelim = 0.44;   // elimination from central
    p.k12 = 0.2;   // central -> peripheral
    p.k21 = 0.1;   // peripheral -> central

    double t_end = 24.0;
    double h = 0.001;

    vector<double> t;
    vector<vector<double>> y;

    three_compartment_model(p, t_end, h, cout, t, y);
    plot_run(t, y, p);

    return 0;
}

// IV multiple bolus dosing (2-compartment)
// Example: every 8 hours.
int run_iv_multidose()
{
    PKParams p;

    p.dose = 10;   // mg each bolus
    p.dose_steps  = 12000; //twice daily
    p.dose_compartment = 0; // central
    p.F = 1.0;

    p.kabsorp = 1.16; //<---------------------------------------------------------------------------------------- Change these
    p.kelim = 0.44;
    p.k12 = 0.2;
    p.k21 = 0.1;

    double t_end = 120;
    double h = 0.001;

    vector<double> t;
    vector<vector<double>> y;

    two_compartment_model(p, t_end, h, cout, t, y);
    plot_run(t, y, p);

    return 0;
}

// Oral multiple dosing (3-compartment)
int run_oral_multidose()
{
    PKParams p;

    p.dose = 10*0.239; // mg each dose
    p.dose_steps = 4000; // twice daily
    p.dose_compartment = 0; // gut
    p.F = 1.;

    p.kabsorp = 1.16; //<---------------------------------------------------------------------------------------- Change these
    p.kelim = 0.44;
    p.k12 = 0.2;
    p.k21 = 0.1;

    double t_end = 48.0;
    double h = 0.001;

    vector<double> t;
    vector<vector<double>> y;

    three_compartment_model(p, t_end, h, cout, t, y);
    plot_run(t, y, p);

    return 0;
}

int run_iv_continuous()
{
    PKParams p;

    p.dose = 0.;
    p.dose_steps = 0;
    p.dose_compartment = 0;
    p.F = 1.0;

    // Continuous infusion rate (mg/h) into central compartment
    p.R_in  = 0.8666;   // e.g. 1 mg/
    p.R_gut = 0.0;

    // Rate constants (1/h)
    p.kabsorp = 1.16; //<---------------------------------------------------------------------------------------- Change these
    p.kelim = 0.44;
    p.k12 = 0.2;
    p.k21 = 0.1;

    double t_end = 120;
    double h = 0.001;

    vector<double> t;
    vector<vector<double>> y;

    two_compartment_model(p, t_end, h, cout, t, y);
    plot_run(t, y, p);

    return 0;
}

int run_oral_continuous()
{
    PKParams p;

    p.dose = 0.0; //zero
    p.dose_steps = 0;
    p.dose_compartment = 0;
    p.F = 1;

    p.R_in  = 0.0;
    p.R_gut = 2.5*0.239; // mg/h zero-order release into gut

    p.kabsorp = 1.16; //<-- need increase
    p.kelim = 0.44;
    p.k12 = 0.2;
    p.k21 = 0.1;

    double t_end = 48.0;
    double h = 0.001;

    vector<double> t;
    vector<vector<double>> y;

    three_compartment_model(p, t_end, h, cout, t, y);
    plot_run(t, y, p);

    return 0;
}



int main()
{
    int choice = 0;

    std::cout << "Select PK simulation mode:\n";
    std::cout << "  1. IV bolus (single)\n";
    std::cout << "  2. IV multidose\n";
    std::cout << "  3. IV continuous infusion\n";
    std::cout << "  4. Oral bolus (single)\n";
    std::cout << "  5. Oral multidose\n";
    std::cout << "  6. Oral continuous release\n";
    std::cout << "  Enter choice (1-6): ";

    std::cin >> choice;
 

    switch (choice)
    {
        case 1:
            return run_iv_bolus();

        case 2:
            return run_iv_multidose();

        case 3:
            return run_iv_continuous();

        case 4:
            return run_oral_bolus();

        case 5:
            return run_oral_multidose();

        case 6:
            return run_oral_continuous();
            
        /*default:
            std::cerr << "Invalid choice.\n";
            return 1;*/
            
        
    }

}
