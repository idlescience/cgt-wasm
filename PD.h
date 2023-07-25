#include "common.h"

void
PD(bool &disp, unsigned short int &n, vector<double> &v, unsigned short int &iter, unsigned int &piv, unsigned int &sr,
   double &t, vector<double> &x, unsigned int &s, bool &nlsu);

void iteration(vector<bool> &unsettled, unsigned int &s, double &xS, unsigned short int &n, vector<vector<bool>> &A,
               vector<double> &x, vector<double> &v, double &epsi, double &prec, vector<vector<double>> &Arref,
               vector<bool> &J, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled,
               vector<double> &settled_values, unsigned short int &iter, unsigned int &sr, vector<bool> &unsettled_p,
               vector<double> &singleton_bounds, bool &nlsu);

void subroutine(vector<bool> &U, vector<bool> &U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2,
                vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n,
                unsigned int &t_size, unsigned short int &t2_size, unsigned short int &rank, bool &disp,
                vector<vector<bool>> &Asettled, unsigned int &sr, vector<double> &settled_values,
                vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi,
                vector<double> &v, vector<unsigned int> &T2_coord);

void
subr_upd(vector<vector<double>> &Arref, vector<bool> &J, unsigned int &i, glp_prob &lp, vector<int> &constr_indices,
         vector<vector<int>> &bal_indices, int (&ia)[], int (&ja)[], double (&ar)[], unsigned short int &n,
         double &prec, vector<bool> &U, vector<bool> &U2, unsigned int &sumt, unsigned short int &sumt2,
         vector<bool> &t, vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2,
         unsigned int &t_size, unsigned short int &t2_size, IloCplex &SR, IloNumVarArray &lambda,
         unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, vector<double> &settled_values,
         vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector<double> &v,
         vector<unsigned int> &T2_coord, IloExpr &sr_obj, IloExprArray &bal_eq, vector<double> &u);
