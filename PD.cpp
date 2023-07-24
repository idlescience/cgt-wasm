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

    // GLPK
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
        const char *col_name = ("x{" + std::to_string(player_index) + "}").c_str();
        glp_set_col_name(lp, player_index, col_name);
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
    int constr_index = 0;

    for (unsigned short int j = 0; j < n; j++) {
        int player_index = j + 1;
        Asettled[0][j] = true;

        // x({i}) >= v({i}); for all i = 1 ... n;
        constr_index++;
        glp_add_rows(lp, 1);
        impu_constr_indices[j] = constr_index;
        const char *row_name = ("x({" + std::to_string(player_index) + "}) >= v({" + std::to_string(j) + "})" + " = " +
                                std::to_string(v[player_index])).c_str();
        glp_set_row_name(lp, constr_index, row_name);
        glp_set_row_bnds(lp, constr_index, GLP_LO, singleton_bounds[j], 0.0);
        ia[ia_index] = constr_index;
        ja[ia_index] = player_index;
        ar[ia_index] = 1.0;
        ia_index++;
    }

    for (unsigned int i = 0; i < s; i++) {
        // x(S) + e(1) >= v(S); S subset of N;
        constr_index++;
        glp_add_rows(lp, 1);
        unsett_ineq_indices[i] = constr_index;
        const char *row_name = "x(S) + e(1) >= v(S)";
        glp_set_row_name(lp, constr_index, row_name);
        glp_set_row_bnds(lp, constr_index, GLP_LO, v[i], 0.0);
        ia[ia_index] = constr_index;
        ja[ia_index] = epsi_index;
        ar[ia_index] = 1.0;
        ia_index++;

        for (unsigned short int j = 0; j < n; j++) {
            int player_index = j + 1;
            if (A[i][j]) {
                ia[ia_index] = constr_index;
                ja[ia_index] = player_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
        }
    }

    // x(N) = v(N)
    constr_index++;
    glp_add_rows(lp, 1);
    eq_indices[s] = constr_index;
    const char *row_name = "x(N) = v(N)";
    glp_set_row_name(lp, constr_index, row_name);
    glp_set_row_bnds(lp, constr_index, GLP_FX, v[s], 0.0);
    for (unsigned short int j = 0; j < n; j++) {
        int player_index = j + 1;
        ia[ia_index] = constr_index;
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
            if (glp_get_row_dual(lp, unsett_ineq_indices[i]) > prec) {
                if (binrank(Arref, J, A[i], n)) {
                    rank++;
                    if (disp) {
                        cout << "Dual: lambda_" << i + 1 << " > 0, rank = " << rank << " (" << s - i
                             << " settled as well)" << endl;
                    }
                    if (rank == n) {
                        t = cpuTime() - t1;
                        if (disp) {
                            cout << "Rank condition satisfied!" << endl;
                        }
                        cout << "finished!" << endl;
                        glp_delete_prob(lp);
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
            if (glp_get_row_dual(lp, impu_constr_indices[j]) > prec) {
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
                        cout << "finished!" << endl;
                        glp_delete_prob(lp);
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
                if (A[i][j]) {
                    xS += x[j];
                }
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
        if (t_size > 0) {
            cout << "Tight coalitions:" << endl;
            for (unsigned short int i = 0; i < t_size; i++) {
                cout << T_coord[i] + 1 << endl;
            }
        }
        if (t2_size > 0) {
            cout << "T0:" << endl;
            for (unsigned short int i = 0; i < t2_size; i++) {
                cout << T2_coord[i] + 1 << endl;
            }
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

    // GLPK
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
        const char *col_name = ("x{" + std::to_string(player_index) + "}").c_str();
        glp_set_col_name(lp, player_index, col_name);
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
    int constr_index = 1;

    for (unsigned short int j = 0; j < n; j++) {
        if (unsettled_p[j]) {
            int player_index = j + 1;

            // x({i}) >= v({i}); for all i = 1 ... n;
            constr_index++;
            glp_add_rows(lp, 1);
            impu_constr_indices[j] = constr_index;
            const char *row_name = ("x({" + std::to_string(player_index) + "}) >= v({" + std::to_string(j) + "})" +
                                    " = " + std::to_string(v[j])).c_str();
            glp_set_row_name(lp, constr_index, row_name);
            glp_set_row_bnds(lp, constr_index, GLP_LO, singleton_bounds[j], 0.0);
            ia[ia_index] = constr_index;
            ja[ia_index] = player_index;
            ar[ia_index] = 1.0;
            ia_index++;
        }
    }

    // x(N) = v(N)
    constr_index++;
    glp_add_rows(lp, 1);
    eq_indices[s] = constr_index;
    const char *row_name = "x(N) = v(N)";
    glp_set_row_name(lp, constr_index, row_name);
    glp_set_row_bnds(lp, constr_index, GLP_FX, v[s], 0.0);
    for (unsigned short int j = 0; j < n; j++) {
        int player_index = j + 1;
        ia[ia_index] = constr_index;
        ja[ia_index] = player_index;
        ar[ia_index] = 1.0;
        ia_index++;
    }

    for (unsigned int i = 0; i < s; i++) {
        if (unsettled[i]) {
            // x(S) + e(k) >= v(S); S in unsettled;
            constr_index++;
            glp_add_rows(lp, 1);
            unsett_ineq_indices[i] = constr_index;
            const char *row_name = "x(S) + e(k) >= v(S)";
            glp_set_row_name(lp, constr_index, row_name);
            glp_set_row_bnds(lp, constr_index, GLP_LO, v[i], 0.0);
            ia[ia_index] = constr_index;
            ja[ia_index] = epsi_index;
            ar[ia_index] = 1.0;
            ia_index++;

            for (unsigned short int j = 0; j < n; j++) {
                int player_index = j + 1;
                if (A[i][j]) {
                    ia[ia_index] = constr_index;
                    ja[ia_index] = player_index;
                    ar[ia_index] = 1.0;
                    ia_index++;
                }
            }
        }
    }

    for (unsigned int i = 1; i < rank; i++) {
        // x(S) = v(S); S in settled;
        constr_index++;
        glp_add_rows(lp, 1);
        unsett_ineq_indices[i] = constr_index;
        const char *row_name = "x(S) + e(1) >= v(S)";
        glp_set_row_name(lp, constr_index, row_name);
        glp_set_row_bnds(lp, constr_index, GLP_FX, settled_values[i], 0.0);
        ia[ia_index] = constr_index;
        ja[ia_index] = epsi_index;
        ar[ia_index] = 0.0;
        ia_index++;

        for (unsigned short int j = 0; j < n; j++) {
            int player_index = j + 1;
            if (Asettled[i][j]) {
                ia[ia_index] = constr_index;
                ja[ia_index] = player_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
        }
    }

    if (disp) {
        cout << endl << "   ---===   SOLVING THE " << iter + 1 << "-TH LP   ===---   " << endl << endl;
    }

    // NO warm start option. Kept going with solving from scratch

    glp_load_matrix(lp, ia_index - 1, ia, ja, ar);
    glp_simplex(lp, NULL);
    iter++;

    for (unsigned short int j = 0; j < n; j++) {
        int player_index = j + 1;
        x[j] = glp_get_col_prim(lp, player_index);
    }

    epsi = glp_get_col_prim(lp, epsi_index);

    if (disp) {
        cout << "New solution point:" << endl;
        for (unsigned short int i = 0; i < n; i++) {
            cout << x[i] << endl;
        }
        cout << "Epsilon: " << epsi << endl;
    }

    for (unsigned int i = 0; i < s; i++) {
        if (unsettled[i]) {
            if (glp_get_row_dual(lp, unsett_ineq_indices[i]) > prec) {
                if (binrank(Arref, J, A[i], n)) {
                    rank++;
                    if (disp) {
                        cout << "Dual: lambda_" << i + 1 << " > 0, rank = " << rank << " (" << s - i
                             << " settled as well)" << endl;
                    }
                    if (rank == n) {
                        if (disp) {
                            cout << "Rank condition satisfied!" << endl;
                        }
                        glp_delete_prob(lp);
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

    for (unsigned short int j = 0; j < n; j++) {
        if (unsettled_p[j]) {
            if (glp_get_row_dual(lp, impu_constr_indices[j]) > prec) {
                if (binrank(Arref, J, A[pow(2, j) - 1], n)) {
                    rank++;
                    if (disp) {
                        cout << "Dual: lambda_impu" << j + 1 << " > 0, rank = " << rank << " (" << s - pow(2, j)
                             << " settled as well)" << endl;
                    }
                    if (rank == n) {
                        if (disp) {
                            cout << "Rank condition satisfied!" << endl;
                        }
                        glp_delete_prob(lp);
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
        cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
    }
}

void subroutine(vector<bool> &U, vector<bool> &U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2,
                vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n,
                unsigned int &t_size, unsigned short int &t2_size, unsigned short int &rank, bool &disp,
                vector<vector<bool>> &Asettled, unsigned int &sr, vector<double> &settled_values,
                vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi,
                vector<double> &v, vector<unsigned int> &T2_coord) {
    // GLPK
    glp_prob *lp;
    int ia[1 + 1000], ja[1 + 1000];
    double ar[1 + 1000];
    lp = glp_create_prob();
    glp_set_prob_name(lp, "SUBROUTINE(k)");
    glp_set_obj_dir(lp, GLP_MAX);

    // objective function
    glp_add_cols(lp, t_size + t2_size);

    for (unsigned short int k = 0; k < t_size + t2_size; k++) {
        // set column names and coefficients
        // sr_obj
        int sr_obj_index = k + 1;
        const char *col_name = ("sr{" + std::to_string(sr_obj_index) + "}").c_str();
        glp_set_col_name(lp, sr_obj_index, col_name);
        glp_set_col_bnds(lp, sr_obj_index, GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(lp, sr_obj_index, 1);
    }

    vector<int> bal_indices(n + 1, -1);
    int ia_index = 1;
    int constr_index = 0;

    for (unsigned short int i = 0; i < n; i++) {
        // eq == 0
        constr_index++;
        glp_add_rows(lp, 1);
        bal_indices[i] = constr_index;
        const char *row_name = ("t({" + std::to_string(i + 1) + "}) = 0").c_str();
        glp_set_row_name(lp, constr_index, row_name);
        glp_set_row_bnds(lp, constr_index, GLP_FX, 0, 0.0);

        for (unsigned int j = 0; j < t_size; j++) {
            if (Atight[j][i] == true) {
                int tight_index = j + 1;
                ia[ia_index] = constr_index;
                ja[ia_index] = tight_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
            if (j < rank && Asettled[j][i] == true) {
                int tight_index = j + t_size + t2_size + 1;
                ia[ia_index] = constr_index;
                ja[ia_index] = tight_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
        }
        for (unsigned int j = 0; j < t2_size; j++) {
            if (Atight2[j][i] == true) {
                int tight_index = j + t_size + 1;
                ia[ia_index] = constr_index;
                ja[ia_index] = tight_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
        }
        if (rank > t_size) {
            for (unsigned short int j = t_size; j < rank; j++) {
                if (Asettled[j][i] == true) {
                    int tight_index = j + t_size + t2_size + 1;
                    ia[ia_index] = constr_index;
                    ja[ia_index] = tight_index;
                    ar[ia_index] = 1.0;
                    ia_index++;
                }
            }
        }
    }

    // pos_eq == 1
    constr_index++;
    glp_add_rows(lp, 1);
    bal_indices[n] = constr_index;
    const char *row_name = ("t({" + std::to_string(n) + "}) = 1").c_str();
    glp_set_row_name(lp, constr_index, row_name);
    glp_set_row_bnds(lp, constr_index, GLP_FX, 1, 0.0);


    vector<int> pos_eq_indices(t_size, -1);

    for (unsigned short int i = 0; i < t_size; i++) {
        int tight_index = i + 1;
        pos_eq_indices[i] = ia_index;
        ia[ia_index] = constr_index;
        ja[ia_index] = tight_index;
        ar[ia_index] = 1.0;
        ia_index++;
    }

    if (disp) {
        cout << endl << "   ---===   SOLVING SUBROUTINE LP   ===---   " << endl << endl;
    }

    glp_load_matrix(lp, ia_index - 1, ia, ja, ar);
    int feas = glp_simplex(lp, NULL);

    if (disp) {
        cout << "subroutine feasibility: " << feas << endl;
    }

    sr++;

    vector<double> u(t_size + t2_size, 0);

    if (feas == 0) {
        for (unsigned short int k = 0; k < t_size + t2_size; k++) {
            const int k_index = k + 1;
            u[k] = glp_get_col_prim(lp, k_index);
        }
    }

    unsigned int i;
    unsigned int sumt = 0;
    unsigned short int sumt2 = 0;
    vector<bool> t(t_size, false);
    vector<bool> t2(t2_size, false);

    while (feas) {
        subr_upd(Arref, J, i, lp, pos_eq_indices, ar, n, prec, U, U2, sumt, sumt2, t, t2, Atight, Atight2, t_size, t2_size, SR,
                 lambda, rank, disp, Asettled, settled_values, unsettled, T_coord, s, epsi, v, T2_coord, sr_obj, bal_eq,
                 u);
        if (rank == n) {
            return;
        } else {
            i = 0;
            while (i < t_size) {
                if (t[i] == false && !(binrank(Arref, J, Atight[i], n))) {
                    U[i] = false;
                    t[i] = true;

                    const int pos_eq_ia_index = pos_eq_indices[i];
                    ar[pos_eq_ia_index] = 0;

                    const int sr_obj_index = i + 1;
                    glp_set_obj_coef(lp, sr_obj_index, 0);

                    for (unsigned short int j = 0; j < n; j++) {
                        if (Atight[i][j]) {
                            bal_eq[j] -= lambda[i];
                        }
                    }
                    sumt++;
                    unsettled[T_coord[i]] = false;
                    unsettled[s - 1 - T_coord[i]] = false;
                    if (disp) {
                        cout << "SETTLED: " << T_coord[i] + 1 << " at " << v[T_coord[i]] - epsi << endl;
                    }
                    if (disp) {
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
                            if (Atight2[i][j]) {
                                bal_eq[j] -= lambda[i + t_size];
                            }
                        }
                        sumt2++;
                        if (disp) {
                            cout << "SETTLED: " << T2_coord[i] + 1 << " at " << v[T2_coord[i]] << endl;
                        }
                        if (unsettled[T2_coord[i]]) {
                            unsettled[T2_coord[i]] = false;
                        }
                        if (unsettled[s - 1 - T2_coord[i]]) {
                            unsettled[s - 1 - T2_coord[i]] = false;
                            if (disp) {
                                cout << T2_coord[i] + 1 << " and " << s - T2_coord[i]
                                     << " got settled without rank increase." << endl;
                            }
                        }
                    }
                }
                i++;
            }
        }
        if (sumt == t_size) {
            return;
        } else {
            sr_model.remove(bal[n]);
            bal_eq[n] = pos_eq;
            r = (pos_eq == 1);
            bal[n] = r;
            sr_model.add(bal[n]);
            sr_model.remove(OBJ);
            OBJ = IloMinimize(sr_env, sr_obj);
            sr_model.add(OBJ);
            if (disp) {
                cout << endl << "   ---===   SOLVING SUBROUTINE LP AGAIN  ===---   " << endl << endl;
            }
            feas = SR.solve();
            if (disp) {
                cout << "subroutine feasibility: " << feas << endl;
            }
            sr++;
            if (feas) {
                for (unsigned short int j = 0; j < t_size + t2_size; j++) {
                    u[j] = SR.getValue(lambda[j]);
                }
            }
        }
    }
    glp_delete_prob(lp);
}

void subr_upd(vector<vector<double>> &Arref, vector<bool> &J, unsigned int &i, int (&ia)[], int (&ja)[], int (&ar)[],
              unsigned short int &n, double &prec, vector<bool> &U, vector<bool> &U2, unsigned int &sumt,
              unsigned short int &sumt2, vector<bool> &t, vector<bool> &t2, vector<vector<bool>> &Atight,
              vector<vector<bool>> &Atight2, unsigned int &t_size, unsigned short int &t2_size, IloCplex &SR,
              IloNumVarArray &lambda, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled,
              vector<double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s,
              double &epsi, vector<double> &v, vector<unsigned int> &T2_coord, IloExpr &sr_obj, IloExprArray &bal_eq,
              vector<double> &u) {

    i = 0;
    while (i < t_size && sumt < t_size) {
        if (t[i] == false && u[i] > prec) {
            U[i] = false;
            t[i] = true;
            pos_eq -= lambda[i];
            sr_obj -= lambda[i];
            sumt++;
            if (binrank(Arref, J, Atight[i], n)) {
                if (disp) {
                    cout << "Rank increased to " << rank + 1 << " with " << T_coord[i] + 1 << " (and " << s - T_coord[i]
                         << ") getting settled." << endl;
                }
                if (rank == n - 1) {
                    rank++;
                    if (disp) {
                        cout << "Rank condition satisfied!" << endl;
                    }
                    return;
                }
                rowechform(Arref, J, Atight[i], n, rank);
                rank++;
                Asettled[rank - 1] = Atight[i];
                settled_values[rank - 1] = v[T_coord[i]] - epsi;
                unsettled[T_coord[i]] = false;
                unsettled[s - 1 - T_coord[i]] = false;
                lambda[i].setLB(-IloInfinity);
                if (disp) {
                    cout << "SETTLED: " << T_coord[i] + 1 << " at " << v[T_coord[i]] - epsi << endl;
                }
            } else {
                for (unsigned short int j = 0; j < n; j++) {
                    if (Atight[i][j]) {
                        bal_eq[j] -= lambda[i];
                    }
                }
                unsettled[T_coord[i]] = false;
                unsettled[s - 1 - T_coord[i]] = false;
                if (disp) {
                    cout << "SETTLED: " << T_coord[i] + 1 << " at " << v[T_coord[i]] - epsi << endl;
                }
                if (disp) {
                    cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase."
                         << endl;
                }
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
                if (disp) {
                    cout << "Rank increased to " << rank + 1 << " with " << T2_coord[i] + 1 << " (and "
                         << s - T2_coord[i] << ") getting settled." << endl;
                }
                if (rank == n - 1) {
                    rank++;
                    if (disp) {
                        cout << "Rank condition satisfied!" << endl;
                    }
                    return;
                }
                rowechform(Arref, J, Atight2[i], n, rank);
                rank++;
                Asettled[rank - 1] == Atight2[i];
                settled_values[rank - 1] = v[T2_coord[i]];
                unsettled[T2_coord[i]] = false;
                unsettled[s - 1 - T2_coord[i]] = false;
                lambda[i + t_size].setLB(-IloInfinity);
                if (disp) {
                    cout << "SETTLED: " << T2_coord[i] + 1 << " at " << v[T2_coord[i]] << endl;
                }
            } else {
                for (unsigned short int j = 0; j < n; j++) {
                    if (Atight2[i][j]) {
                        bal_eq[j] -= lambda[i + t_size];
                    }
                }
                unsettled[T2_coord[i]] = false;
                unsettled[s - 1 - T2_coord[i]] = false;
                if (disp) {
                    cout << "SETTLED: " << T2_coord[i] + 1 << " at " << v[T2_coord[i]] << endl;
                }
                if (disp) {
                    cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase."
                         << endl;
                }
            }
        }
        i++;
    }
}