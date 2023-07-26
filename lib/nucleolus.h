#include "common.h"

extern "C" void PD(double *v_orig, double *x_orig, unsigned short int n);

void iteration(vector<bool> &unsettled, unsigned int &s, double &xS, unsigned short int &n, vector<vector<bool>> &A,
               double *x, double *v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J,
               int &rank, bool &disp, vector<vector<bool>> &Asettled, vector<double> &settled_values,
               unsigned short int &iter, unsigned int &sr, vector<bool> &unsettled_p, vector<double> &singleton_bounds,
               bool &nlsu);

void subroutine(vector<bool> &U, vector<bool> &U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2,
                vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n, int &t_size,
                int &t2_size, int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sr,
                vector<double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s,
                double &epsi, double *v, vector<unsigned int> &T2_coord);

void subr_upd(vector<vector<double>> &Arref, vector<bool> &J, glp_prob *lp, int &ia_index, vector<int> constr_indices,
              vector<vector<int>> bal_indices, int ia[], int ja[], double ar[], unsigned short int &n, double &prec,
              vector<bool> &U, vector<bool> &U2, unsigned int &sumt, unsigned short int &sumt2, vector<bool> &t,
              vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, int &t_size, int &t2_size,
              int &rank, bool &disp, vector<vector<bool>> &Asettled, vector<double> &settled_values,
              vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, double *v,
              vector<unsigned int> &T2_coord, vector<double> &u);