/*
 *    Nucleolus
 *    PD.cpp
 *    Purpose: finding the nucleolus-wasm of a cooperative game using the
 *             primal-dual sequence of Benedek (2019) - Computing the
 *             nucleolus-wasm of cooperative games
 *
 *    @author Marton Benedek
 *    @version 1.0 16/07/2019
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program. If not, see:
 *    <https://github.com/blrzsvrzs/nucleolus>.
 */

#include "PD.h"

void PD(bool &disp, unsigned short int &n, vector<double> &v, unsigned short int &iter, unsigned int &sr, double &t,
        vector<double> &x, unsigned int &s, bool &nlsu) {
    double prec = pow(10, -6);
    vector<bool> unsettled(s + 1, true);
    unsettled[s] = false;
    double t1 = cpuTime();
    vector<vector<bool>> A(s + 1, vector<bool>(n, false));
    A_mx(A, n, s);
    vector<double> singleton_bounds(n, 0);
    double impu = 0;

    for (unsigned short int i = 0; i < n; i++) {
        singleton_bounds[i] = v[pow(2, i) - 1];
        impu += singleton_bounds[i];
    }

    for (unsigned short int i = 0; i < n; i++) {
        x[i] = singleton_bounds[i] + (v[s] - impu) / n;
    }

    vector<bool> unsettled_p(n, true);
    vector<vector<double>> Arref(n, vector<double>(n, 0));
    Arref[0] = vector<double>(n, 1);
    vector<bool> J(n, true);
    J[0] = false;
    unsigned short int rank = 1;
    vector<vector<bool>> Asettled(n, vector<bool>(n, false));
    vector<double> settled_values(n, 0);
    settled_values[0] = v[s];
    double epsi = 0;

    // IloEnv env;
    glp_prob *lp;
    int ia[1 + 1000], ja[1 + 1000];
    double ar[1 + 1000];
    lp = glp_create_prob();
    glp_set_prob_name(lp, "P(1)");
    glp_set_obj_dir(lp, GLP_MIN);

    // objective function
    glp_add_cols(lp, n + 1);

    int ia_index = 1;

    for (unsigned short int j = 0; j < n; j++) {
        // set column names and coefficients
        int player_index = j + 1;
        const char *player_col_name = ("x{" + std::to_string(j) + "}").c_str();
        glp_set_col_name(lp, player_index, player_col_name);
        glp_set_col_bnds(lp, player_index, GLP_FR, 0.0, 0.0);
        glp_set_obj_coef(lp, player_index, 0.0);
    }

    // the only coefficient of the objective function, epsi
    int epsi_index = n + 1;
    glp_set_col_name(lp, epsi_index, "e(1)");
    glp_set_col_bnds(lp, epsi_index, GLP_FR, 0.0, 0.0);
    glp_set_obj_coef(lp, epsi_index, 1);

    vector<int> eq_indices(s, -1);
    vector<int> unsett_ineq_indices(s, -1);
    vector<int> impu_constr_indices(n, -1);
    int constraintIndex = 0;

    for (unsigned short int j = 0; j < n; j++) {
        int player_index = j + 1;
        Asettled[0][j] = true;

        // x({i}) >= v({i}); for all i = 1 ... n;
        constraintIndex++;
        glp_add_rows(lp, 1);
        impu_constr_indices[j] = constraintIndex;
        const char *individual_constraint_row_name = ("x({" + std::to_string(j) + "}) >= v({" + std::to_string(j) +
                                                      "})" + " = " + std::to_string(v[j])).c_str();
        glp_set_row_name(lp, constraintIndex, individual_constraint_row_name);
        glp_set_row_bnds(lp, constraintIndex, GLP_LO, singleton_bounds[j], 0.0);
        ia[ia_index] = constraintIndex;
        ja[ia_index] = player_index;
        ar[ia_index] = 1.0;
        ia_index++;
    }

    for (unsigned int i = 0; i < s; i++) {
        // x(S) + e(1) >= v(S); S subset of N;
        constraintIndex++;
        glp_add_rows(lp, 1);
        unsett_ineq_indices[i] = constraintIndex;
        const char *excess_constraint_col_name = "x(S) + e(1) >= v(S)";
        glp_set_row_name(lp, constraintIndex, excess_constraint_col_name);
        glp_set_row_bnds(lp, constraintIndex, GLP_LO, v[i], 0.0);
        ia[ia_index] = constraintIndex;
        ja[ia_index] = epsi_index;
        ar[ia_index] = 1.0;
        ia_index++;

        for (unsigned short int j = 0; j < n; j++) {
            int player_index = j + 1;
            if (A[i][j]) {
                ia[ia_index] = constraintIndex;
                ja[ia_index] = player_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
        }
    }

    // x(N) = v(N)
    constraintIndex++;
    glp_add_rows(lp, 1);
    eq_indices[s] = constraintIndex;
    const char *eq_constraint_row_name = "x(N) = v(N)";
    glp_set_row_name(lp, constraintIndex, eq_constraint_row_name);
    glp_set_row_bnds(lp, constraintIndex, GLP_FX, v[s], 0.0);
    for (unsigned short int j = 0; j < n; j++) {
        int player_index = j + 1;
        ia[ia_index] = constraintIndex;
        ja[ia_index] = player_index;
        ar[ia_index] = 1.0;
        ia_index++;
    }

    glp_load_matrix(lp, ia_index - 1, ia, ja, ar);
    glp_simplex(lp, NULL);
    iter++;

    for (unsigned short int j = 0; j < n; j++) {
        int player_index = j + 1;
        x[j] = glp_get_col_prim(lp, player_index);
    }

    epsi = glp_get_col_prim(lp, epsi_index);
    if (disp) {
        cout << "Least core solution:" << endl;
        for (unsigned short int j = 0; j < n; j++) {
            cout << x[j] << endl;
        }
        cout << "Least core value: " << epsi << endl;
    }
    for (unsigned int i = 0; i < s; i++) {
        if (unsettled[i]) {
            int coalition_index = unsett_ineq_indices[i];
            if (glp_get_row_dual(lp, coalition_index) > prec) {
                if (binrank(Arref, J, A[i], n)) {
                    rank++;
                    if (disp)
                        cout << "Dual: lambda_" << i + 1 << " > 0, rank = " << rank << " (" << s - i
                             << " settled as well)" << endl;
                    if (rank == n) {
                        t = cpuTime() - t1;
                        if (disp) {
                            cout << "Rank condition satisfied!" << endl;
                        }
                        glp_delete_prob(lp);
                        cout << "finished!" << endl;
                        return;
                    }
                    rowechform(Arref, J, A[i], n, rank);
                    Asettled[rank - 1] = A[i];
                    settled_values[rank - 1] = v[i] - epsi;
                    if (disp)
                        cout << "SETTLED: " << i + 1 << " at " << v[i] - epsi << endl;
                }
                unsettled[i] = false;
                unsettled[s - 1 - i] = false;
            }
        }
    }

    for (unsigned short int j = 0; j < n; j++) {
        if (unsettled_p[j]) {
            int player_index = impu_constr_indices[j];
            if (glp_get_row_dual(lp, player_index) > prec) {
                if (binrank(Arref, J, A[pow(2, j) - 1], n)) {
                    rank++;
                    if (disp) {
                        cout << "Dual: lambda_impu" << j + 1 << " > 0, rank = " << rank << " (" << s - pow(2, j)
                             << " settled as well)" << endl;
                    }
                    if (rank == n) {
                        t = cpuTime() - t1;
                        if (disp) {
                            cout << "Rank condition satisfied!" << endl;
                        }
                        glp_delete_prob(lp);
                        cout << "finished!" << endl;
                        return;
                    }
                    rowechform(Arref, J, A[pow(2, j) - 1], n, rank);
                    Asettled[rank - 1] = A[pow(2, j) - 1];
                    settled_values[rank - 1] = v[pow(2, j) - 1];
                    if (disp) {
                        cout << "SETTLED: " << pow(2, j) << " at " << v[pow(2, j) - 1] << endl;
                    }
                }
                unsettled[pow(2, j) - 1] = false;
                unsettled[s - pow(2, j)] = false;
                unsettled_p[j] = false;
            }
        }
    }

    glp_delete_prob(lp);

    if (disp) {
        cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
    }
    double xS;
    while (rank < n) {
        iteration(unsettled, s, xS, n, A, x, v, epsi, prec, Arref, J, rank, disp, Asettled, settled_values, iter, sr,
                  unsettled_p, singleton_bounds, nlsu);
    }
    t = cpuTime() - t1;
    cout << "PD finished!" << endl;
    cout << "The nucleolus-wasm solution:" << endl;
    for (unsigned short int i = 0; i < n; i++)
        cout << x[i] << endl;
    cout << "Time needed: " << t << " seconds" << endl;
    cout << "Iterations needed: " << iter << endl;
    cout << "Subroutine solves needed: " << sr << endl;
}

void iteration(vector<bool> &unsettled, unsigned int &s, double &xS, unsigned short int &n, vector<vector<bool>> &A,
               vector<double> &x, vector<double> &v, double &epsi, double &prec, vector<vector<double>> &Arref,
               vector<bool> &J, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled,
               vector<double> &settled_values, unsigned short int &iter, unsigned int &sr, vector<bool> &unsettled_p,
               vector<double> &singleton_bounds, bool &nlsu) {

    vector<bool> T(s, false);
    vector<unsigned int> T_coord(0, 0);
    unsigned int t_size = 0;
    vector<bool> T2(n, false);
    vector<unsigned int> T2_coord(0, 0);
    unsigned short int t2_size = 0;

    for (unsigned int i = 0; i < s; i++) {
        if (unsettled[i]) {
            xS = 0;
            for (unsigned short int j = 0; j < n; j++) {
                if (A[i][j])
                    xS += x[j];
            }
            if (abs(v[i] - xS - epsi) < prec) {
                T[i] = true;
                T_coord.push_back(i);
                t_size++;
            }
        }
    }

    for (unsigned int i = 0; i < n; i++) {
        if (unsettled_p[i]) {
            if (abs(x[i] - singleton_bounds[i]) < prec) {
                T2[i] = true;
                T2_coord.push_back(pow(2, i) - 1);
                t2_size++;
            }
        }
    }

    if (disp) {
        cout << "Tight coalitions:" << endl;
        for (unsigned short int i = 0; i < t_size; i++)
            cout << T_coord[i] + 1 << endl;
        if (t2_size > 0) {
            cout << "T0:" << endl;
            for (unsigned short int i = 0; i < t_size; i++)
                cout << T2_coord[i] + 1 << endl;
        }
    }

    vector<vector<bool>> Atight(t_size, vector<bool>(n, false));
    for (unsigned int i = 0; i < t_size; i++) {
        Atight[i] = A[T_coord[i]];
    }

    vector<vector<bool>> Atight2(t2_size, vector<bool>(n, false));
    for (unsigned int i = 0; i < t2_size; i++) {
        Atight2[i] = A[T2_coord[i]];
    }

    vector<bool> U(t_size, true);
    vector<bool> U2(t2_size, true);
    subroutine(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, Asettled, sr, settled_values,
               unsettled, T_coord, s, epsi, v, T2_coord);
    if (disp) {
        cout << endl << "   ---===   SUBROUTINE FINISHED   ===---   " << endl << endl;
        cout << "MIN TIGHT SET FOUND!" << endl;
        for (unsigned int i = 0; i < t_size; i++) {
            if (!U[i])
                cout << T_coord[i] + 1 << endl;
        }
        cout << endl;
    }
    if (rank == n) {
        return;
    }
    if (disp) {
        cout << "Rank increased to: " << rank << endl;
    }
    if (!nlsu) {
        for (unsigned int i = 0; i < s; i++) {
            if (unsettled[i]) {
                if (!(binrank(Arref, J, A[i], n))) {
                    unsettled[i] = false;
                    unsettled[s - 1 - i] = false;
                }
            }
        }
    }

    for (unsigned short int i = 0; i < n; i++) {
        if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false) {
            unsettled_p[i] = false;
        }
    }

    // IloEnv env;
    glp_prob *lp;
    int ia[1 + 1000], ja[1 + 1000];
    double ar[1 + 1000];
    lp = glp_create_prob();
    glp_set_prob_name(lp, "P(2)");
    glp_set_obj_dir(lp, GLP_MIN);

    // objective function
    glp_add_cols(lp, n + 1);

    int ia_index = 1;

    for (unsigned short int j = 0; j < n; j++) {
        // set column names and coefficients
        int player_index = j + 1;
        const char *player_col_name = ("x{" + std::to_string(j) + "}").c_str();
        glp_set_col_name(lp, player_index, player_col_name);
        glp_set_col_bnds(lp, player_index, GLP_FR, 0.0, 0.0);
        glp_set_obj_coef(lp, player_index, 0.0);
    }

    // the only coefficient of the objective function, epsi
    int epsi_index = n + 1;
    glp_set_col_name(lp, epsi_index, "e(1)");
    glp_set_col_bnds(lp, epsi_index, GLP_FR, 0.0, 0.0);
    glp_set_obj_coef(lp, epsi_index, 1);

    vector<int> eq_indices(s, -1);
    vector<int> unsett_ineq_indices(s, -1);
    vector<int> impu_constr_indices(n, -1);
    int constraintIndex = 1;

    for (unsigned short int j = 0; j < n; j++) {
        if (unsettled_p[j]) {
            int player_index = j + 1;

            // x({i}) >= v({i}); for all i = 1 ... n;
            constraintIndex++;
            glp_add_rows(lp, 1);
            impu_constr_indices[j] = constraintIndex;
            const char *individual_constraint_row_name = ("x({" + std::to_string(j) + "}) >= v({" + std::to_string(j) +
                                                          "})" + " = " + std::to_string(v[j])).c_str();
            glp_set_row_name(lp, constraintIndex, individual_constraint_row_name);
            glp_set_row_bnds(lp, constraintIndex, GLP_LO, singleton_bounds[j], 0.0);
            ia[ia_index] = constraintIndex;
            ja[ia_index] = player_index;
            ar[ia_index] = 1.0;
            ia_index++;
        }
    }

    // x(N) = v(N)
    constraintIndex++;
    glp_add_rows(lp, 1);
    eq_indices[s] = constraintIndex;
    const char *eq_constraint_row_name = "x(N) = v(N)";
    glp_set_row_name(lp, constraintIndex, eq_constraint_row_name);
    glp_set_row_bnds(lp, constraintIndex, GLP_FX, v[s], 0.0);
    for (unsigned short int j = 0; j < n; j++) {
        int player_index = j + 1;
        ia[ia_index] = constraintIndex;
        ja[ia_index] = player_index;
        ar[ia_index] = 1.0;
        ia_index++;
    }

    for (unsigned int i = 0; i < s; i++) {
        if (unsettled[i]) {
            // x(S) + e(1) >= v(S); S in unsettled;
            constraintIndex++;
            glp_add_rows(lp, 1);
            unsett_ineq_indices[i] = constraintIndex;
            const char *excess_constraint_col_name = "x(S) + e(1) >= v(S)";
            glp_set_row_name(lp, constraintIndex, excess_constraint_col_name);
            glp_set_row_bnds(lp, constraintIndex, GLP_LO, v[i], 0.0);
            ia[ia_index] = constraintIndex;
            ja[ia_index] = epsi_index;
            ar[ia_index] = 1.0;
            ia_index++;

            for (unsigned short int j = 0; j < n; j++) {
                int player_index = j + 1;
                if (A[i][j]) {
                    ia[ia_index] = constraintIndex;
                    ja[ia_index] = player_index;
                    ar[ia_index] = 1.0;
                    ia_index++;
                }
            }
        }
    }

    for (unsigned int i = 1; i < rank; i++) {
        // x(S) = v(S); S in settled;
        constraintIndex++;
        glp_add_rows(lp, 1);
        unsett_ineq_indices[i] = constraintIndex;
        const char *excess_constraint_col_name = "x(S) + e(1) >= v(S)";
        glp_set_row_name(lp, constraintIndex, excess_constraint_col_name);
        glp_set_row_bnds(lp, constraintIndex, GLP_FX, settled_values[i], 0.0);
        ia[ia_index] = constraintIndex;
        ja[ia_index] = epsi_index;
        ar[ia_index] = 0.0;
        ia_index++;

        for (unsigned short int j = 0; j < n; j++) {
            int player_index = j + 1;
            if (Asettled[i][j]) {
                ia[ia_index] = constraintIndex;
                ja[ia_index] = player_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
        }
    }

    glp_load_matrix(lp, ia_index - 1, ia, ja, ar);
    glp_simplex(lp, NULL);

    if (disp) {
        cout << endl << "   ---===   SOLVING THE " << iter + 1 << "-TH LP   ===---   " << endl << endl;
    }

    IloNumArray x_val(env, n + 1);
    for (unsigned short int i = 0; i < n; i++) {
        x_val[i] = x[i];
    }
    x_val[n] = epsi;
    lp.setStart(x_val, NULL, X, NULL, NULL, NULL);
    lp.solve();
    iter++;
    for (unsigned short int i = 0; i < n; i++) {
        x[i] = lp.getValue(X[i]);
    }
    epsi = lp.getValue(X[n]);
    if (disp) {
        cout << "New solution point:" << endl;
        for (unsigned short int i = 0; i < n; i++) {
            cout << x[i] << endl;
        }
        cout << "Epsilon: " << epsi << endl;
    }

    for (unsigned int i = 0; i < s; i++) {
        if (unsettled[i]) {
            if (lp.getDual(unsett_ineq[i]) > prec) {
                if (binrank(Arref, J, A[i], n)) {
                    rank++;
                    if (disp) {
                        cout << "Dual: lambda_" << i + 1 << " > 0, rank = " << rank << " (" << s - i
                             << " settled as well)" << endl;
                    }
                    if (rank == n) {
                        env.end();
                        if (disp) {
                            cout << "Rank condition satisfied!" << endl;
                        }
                        return;
                    }
                    rowechform(Arref, J, A[i], n, rank);
                    Asettled[rank - 1] = A[i];
                    settled_values[rank - 1] = v[i] - epsi;
                    if (disp) {
                        cout << "SETTLED: " << i + 1 << " at " << v[i] - epsi << endl;
                    }
                }
                unsettled[i] = false;
                unsettled[s - 1 - i] = false;
            }
        }
    }

    for (unsigned short int i = 0; i < n; i++) {
        if (unsettled_p[i]) {
            if (lp.getDual(impu_constr[i]) > prec) {
                if (binrank(Arref, J, A[pow(2, i) - 1], n)) {
                    rank++;
                    if (disp) {
                        cout << "Dual: lambda_impu" << i + 1 << " > 0, rank = " << rank << " (" << s - pow(2, i)
                             << " settled as well)" << endl;
                    }
                    if (rank == n) {
                        env.end();
                        if (disp) {
                            cout << "Rank condition satisfied!" << endl;
                        }
                        return;
                    }
                    rowechform(Arref, J, A[pow(2, i) - 1], n, rank);
                    Asettled[rank - 1] = A[pow(2, i) - 1];
                    settled_values[rank - 1] = v[pow(2, i) - 1];
                    if (disp) {
                        cout << "SETTLED: " << pow(2, i) << " at " << v[pow(2, i) - 1] << endl;
                    }
                }
                unsettled[pow(2, i) - 1] = false;
                unsettled[s - pow(2, i)] = false;
                unsettled_p[i] = false;
            }
        }
    }

    env.end();
    if (disp) {
        cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
    }
}

void subroutine(vector<bool> &U, vector<bool> &U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2,
                vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n,
                unsigned int &t_size, unsigned short int &t2_size, unsigned short int &rank, bool &disp,
                vector<vector<bool>> &Asettled, unsigned int &sr, vector<double> &settled_values,
                vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi,
                vector<double> &v, vector<unsigned int> &T2_coord) {
    unsigned int sumt = 0;
    vector<bool> t(t_size, false);
    unsigned short int sumt2 = 0;
    vector<bool> t2(t2_size, false);

    IloEnv sr_env;
    IloModel sr_model(sr_env);
    IloNumVarArray lambda(sr_env, t_size + t2_size + rank, 0, IloInfinity);
    vector<double> u(t_size + t2_size, 0);
    IloExpr sr_obj(sr_env);
    IloRangeArray bal(sr_env, n + 1);
    IloExprArray bal_eq(sr_env, n + 1);
    for (unsigned int i = 0; i < rank; i++)
        lambda[i + t_size + t2_size].setLB(-IloInfinity);
    IloExpr pos_eq(sr_env);
    for (unsigned short int i = 0; i < n; i++) {
        IloExpr eq(sr_env);
        for (unsigned int j = 0; j < t_size; j++) {
            if (Atight[j][i] == true)
                eq += lambda[j];
            if (i == 0) {
                pos_eq += lambda[j];
                sr_obj += lambda[j];
            }
            if (j < rank && Asettled[j][i] == true)
                eq += lambda[j + t_size + t2_size];
        }
        for (unsigned int j = 0; j < t2_size; j++) {
            if (Atight2[j][i] == true)
                eq += lambda[j + t_size];
            if (i == 0)
                sr_obj += lambda[j + t_size];
        }
        if (rank > t_size) {
            for (unsigned short int j = t_size; j < rank; j++) {
                if (Asettled[j][i] == true)
                    eq += lambda[j + t_size + t2_size];
            }
        }
        bal_eq[i] = eq;
        IloRange r = (eq == 0);
        bal[i] = r;
        sr_model.add(bal[i]);
    }
    bal_eq[n] = pos_eq;
    IloRange r = (pos_eq == 1);
    bal[n] = r;
    sr_model.add(bal[n]);
    IloObjective OBJ = IloMaximize(sr_env, sr_obj);
    sr_model.add(OBJ);
    IloCplex SR(sr_model);
    SR.setParam(IloCplex::Param::RootAlgorithm, 1);
    if (!disp)
        SR.setOut(sr_env.getNullStream());
    else
        cout << endl << "   ---===   SOLVING SUBROUTINE LP   ===---   " << endl << endl;
    bool feas = SR.solve();
    if (disp)
        cout << "subroutine feasibility: " << feas << endl;
    sr++;
    if (feas) {
        for (unsigned short int j = 0; j < t_size + t2_size; j++)
            u[j] = SR.getValue(lambda[j]);
    }
    unsigned int i;
    while (feas) {
        subr_upd(Arref, J, i, pos_eq, n, prec, U, U2, sumt, sumt2, t, t2, Atight, Atight2, t_size, t2_size, SR, lambda,
                 rank, disp, Asettled, settled_values, unsettled, T_coord, s, epsi, v, T2_coord, sr_obj, bal_eq, u);
        if (rank == n)
            return;
        else {
            i = 0;
            while (i < t_size) {
                if (t[i] == false) {
                    if (!(binrank(Arref, J, Atight[i], n))) {
                        U[i] = false;
                        t[i] = true;
                        pos_eq -= lambda[i];
                        sr_obj -= lambda[i];
                        for (unsigned short int j = 0; j < n; j++) {
                            if (Atight[i][j])
                                bal_eq[j] -= lambda[i];
                        }
                        sumt++;
                        unsettled[T_coord[i]] = false;
                        unsettled[s - 1 - T_coord[i]] = false;
                        if (disp)
                            cout << "SETTLED: " << T_coord[i] + 1 << " at " << v[T_coord[i]] - epsi << endl;
                        if (disp)
                            cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase."
                                 << endl;
                    }
                }
                i++;
            }
            i = 0;
            while (i < t2_size) {
                if (t2[i] == false) {
                    if (!(binrank(Arref, J, Atight2[i], n))) {
                        U2[i] = false;
                        t2[i] = true;
                        sr_obj -= lambda[i + t_size];
                        for (unsigned short int j = 0; j < n; j++) {
                            if (Atight2[i][j])
                                bal_eq[j] -= lambda[i + t_size];
                        }
                        sumt2++;
                        if (disp)
                            cout << "SETTLED: " << T2_coord[i] + 1 << " at " << v[T2_coord[i]] << endl;
                        if (unsettled[T2_coord[i]]) {
                            unsettled[T2_coord[i]] = false;
                        }
                        if (unsettled[s - 1 - T2_coord[i]]) {
                            unsettled[s - 1 - T2_coord[i]] = false;
                            if (disp)
                                cout << T2_coord[i] + 1 << " and " << s - T2_coord[i]
                                     << " got settled without rank increase." << endl;
                        }
                    }
                }
                i++;
            }
        }
        if (sumt == t_size)
            return;
        else {
            sr_model.remove(bal[n]);
            bal_eq[n] = pos_eq;
            r = (pos_eq == 1);
            bal[n] = r;
            sr_model.add(bal[n]);
            sr_model.remove(OBJ);
            OBJ = IloMinimize(sr_env, sr_obj);
            sr_model.add(OBJ);
            if (disp)
                cout << endl << "   ---===   SOLVING SUBROUTINE LP AGAIN  ===---   " << endl << endl;
            feas = SR.solve();
            if (disp)
                cout << "subroutine feasibility: " << feas << endl;
            sr++;
            if (feas) {
                for (unsigned short int j = 0; j < t_size + t2_size; j++)
                    u[j] = SR.getValue(lambda[j]);
            }
        }
    }
    sr_env.end();
}

void subr_upd(vector<vector<double>> &Arref, vector<bool> &J, unsigned int &i, IloExpr &pos_eq, unsigned short int &n,
              double &prec, vector<bool> &U, vector<bool> &U2, unsigned int &sumt, unsigned short int &sumt2,
              vector<bool> &t, vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2,
              unsigned int &t_size, unsigned short int &t2_size, IloCplex &SR, IloNumVarArray &lambda,
              unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, vector<double> &settled_values,
              vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector<double> &v,
              vector<unsigned int> &T2_coord, IloExpr &sr_obj, IloExprArray &bal_eq, vector<double> &u) {

    i = 0;
    while (i < t_size && sumt < t_size) {
        if (t[i] == false && u[i] > prec) {
            U[i] = false;
            t[i] = true;
            pos_eq -= lambda[i];
            sr_obj -= lambda[i];
            sumt++;
            if (binrank(Arref, J, Atight[i], n)) {
                if (disp)
                    cout << "Rank increased to " << rank + 1 << " with " << T_coord[i] + 1 << " (and " << s - T_coord[i]
                         << ") getting settled." << endl;
                if (rank == n - 1) {
                    rank++;
                    if (disp)
                        cout << "Rank condition satisfied!" << endl;
                    return;
                }
                rowechform(Arref, J, Atight[i], n, rank);
                rank++;
                Asettled[rank - 1] = Atight[i];
                settled_values[rank - 1] = v[T_coord[i]] - epsi;
                unsettled[T_coord[i]] = false;
                unsettled[s - 1 - T_coord[i]] = false;
                lambda[i].setLB(-IloInfinity);
                if (disp)
                    cout << "SETTLED: " << T_coord[i] + 1 << " at " << v[T_coord[i]] - epsi << endl;
            } else {
                for (unsigned short int j = 0; j < n; j++) {
                    if (Atight[i][j])
                        bal_eq[j] -= lambda[i];
                }
                unsettled[T_coord[i]] = false;
                unsettled[s - 1 - T_coord[i]] = false;
                if (disp)
                    cout << "SETTLED: " << T_coord[i] + 1 << " at " << v[T_coord[i]] - epsi << endl;
                if (disp)
                    cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase."
                         << endl;
            }
        }
        i++;
    }
    i = 0;
    while (i < t2_size && sumt2 < t2_size) {
        if (t2[i] == false && u[i + t_size] > prec) {
            U2[i] = false;
            t2[i] = true;
            sumt2++;
            sr_obj -= lambda[i + t_size];
            if (binrank(Arref, J, Atight2[i], n)) {
                if (disp)
                    cout << "Rank increased to " << rank + 1 << " with " << T2_coord[i] + 1 << " (and "
                         << s - T2_coord[i] << ") getting settled." << endl;
                if (rank == n - 1) {
                    rank++;
                    if (disp)
                        cout << "Rank condition satisfied!" << endl;
                    return;
                }
                rowechform(Arref, J, Atight2[i], n, rank);
                rank++;
                Asettled[rank - 1] == Atight2[i];
                settled_values[rank - 1] = v[T2_coord[i]];
                unsettled[T2_coord[i]] = false;
                unsettled[s - 1 - T2_coord[i]] = false;
                lambda[i + t_size].setLB(-IloInfinity);
                if (disp)
                    cout << "SETTLED: " << T2_coord[i] + 1 << " at " << v[T2_coord[i]] << endl;
            } else {
                for (unsigned short int j = 0; j < n; j++) {
                    if (Atight2[i][j])
                        bal_eq[j] -= lambda[i + t_size];
                }
                unsettled[T2_coord[i]] = false;
                unsettled[s - 1 - T2_coord[i]] = false;
                if (disp)
                    cout << "SETTLED: " << T2_coord[i] + 1 << " at " << v[T2_coord[i]] << endl;
                if (disp)
                    cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase."
                         << endl;
            }
        }
        i++;
    }
}