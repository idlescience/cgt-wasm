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

#include "nucleolus.h"

extern "C" EMSCRIPTEN_KEEPALIVE void PD(double * v_orig, double * x_orig, unsigned short int  n) {
    unsigned int s = (int)(pow(2, n) - 2);
    std::vector<double> v(v_orig, v_orig + s + 1);
    vector<double> x(x_orig, x_orig + n);
    unsigned short int iter = 0;
    unsigned int sr = 0;
    double t = 0;
    bool nlsu = false;
    double prec = pow(10, -6);
    vector<bool> unsettled(s + 1, true);
    unsettled[s] = false;
    double t1 = cpuTime();
    vector<vector<bool>> A(s + 1, vector<bool>(n, false));
    A_mx(A, n, s);
    vector<double> singleton_bounds(n, 0);
    double impu = 0;

    for (int i = 0; i < n; i++) {
        singleton_bounds[i] = v[(long) (pow(2, (double) i) - 1)];
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
    int rank = 1;
    vector<vector<bool>> Asettled(n, vector<bool>(n, false));
    vector<double> settled_values(n, 0);
    settled_values[0] = v[s];
    double epsi = 0;

    // GLPK
    glp_prob *lp;
    int ia[1 + 1000] = {0};
    int ja[1 + 1000] = {0};
    double ar[1 + 1000] = {0};
    lp = glp_create_prob();
    glp_set_prob_name(lp, "P(1)");
    glp_set_obj_dir(lp, GLP_MIN);

    // objective function
    glp_add_cols(lp, n + 1);

    int ia_index = 1;

    for (unsigned short int j = 0; j < n; j++) {
        // set column names and coefficients
        const int col_index = j + 1;
        const std::string col_name_string = "x{" + std::to_string(col_index) + "}";
        const char *col_name = col_name_string.c_str();
        glp_set_col_name(lp, col_index, col_name);
        glp_set_col_bnds(lp, col_index, GLP_FR, 0.0, 0.0);
        glp_set_obj_coef(lp, col_index, 0.0);
    }

    // the only coefficient of the objective function, epsi
    int epsi_index = n + 1;
    glp_set_col_name(lp, epsi_index, "e(1)");
    glp_set_col_bnds(lp, epsi_index, GLP_FR, 0.0, 0.0);
    glp_set_obj_coef(lp, epsi_index, 1);

    vector<int> unsett_ineq_indices(s, -1);
    vector<int> impu_constr_indices(n, -1);
    int constr_index = 0;

    for (unsigned short int j = 0; j < n; j++) {
        const int col_index = j + 1;
        Asettled[0][j] = true;

        // x({i}) >= v({i}); for all i = 1 ... n;
        constr_index++;
        glp_add_rows(lp, 1);
        impu_constr_indices[j] = constr_index;
        const std::string row_name_string =
                "x({" + std::to_string(col_index) + "}) >= v({" + std::to_string(j) + "})" + " = " +
                std::to_string(v[col_index]);
        const char *row_name = row_name_string.c_str();
        glp_set_row_name(lp, constr_index, row_name);
        glp_set_row_bnds(lp, constr_index, GLP_LO, singleton_bounds[j], 0.0);
        ia[ia_index] = constr_index;
        ja[ia_index] = col_index;
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
            const int col_index = j + 1;
            if (A[i][j]) {
                ia[ia_index] = constr_index;
                ja[ia_index] = col_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
        }
    }

    // x(N) = v(N)
    constr_index++;
    glp_add_rows(lp, 1);
    const char *row_name = "x(N) = v(N)";
    glp_set_row_name(lp, constr_index, row_name);
    glp_set_row_bnds(lp, constr_index, GLP_FX, v[s], 0.0);
    for (unsigned short int j = 0; j < n; j++) {
        const int col_index = j + 1;
        ia[ia_index] = constr_index;
        ja[ia_index] = col_index;
        ar[ia_index] = 1.0;
        ia_index++;
    }

    glp_load_matrix(lp, ia_index - 1, ia, ja, ar);
    glp_simplex(lp, nullptr);
    iter++;

    for (unsigned short int j = 0; j < n; j++) {
        const int col_index = j + 1;
        x[j] = glp_get_col_prim(lp, col_index);
    }

    epsi = glp_get_col_prim(lp, epsi_index);
    for (unsigned int i = 0; i < s; i++) {
        if (unsettled[i]) {
            if (glp_get_row_dual(lp, unsett_ineq_indices[i]) > prec) {
                if (binrank(Arref, J, A[i], n)) {
                    rank++;
                    if (rank == n) {
                        t = cpuTime() - t1;
                        cout << "finished!" << endl;
                        glp_delete_prob(lp);
                        return;
                    }
                    rowechform(Arref, J, A[i], n, rank);
                    Asettled[rank - 1] = A[i];
                    settled_values[rank - 1] = v[i] - epsi;
                }
                unsettled[i] = false;
                unsettled[s - 1 - i] = false;
            }
        }
    }

    for (unsigned short int j = 0; j < n; j++) {
        if (unsettled_p[j]) {
            if (glp_get_row_dual(lp, impu_constr_indices[j]) > prec) {
                if (binrank(Arref, J, A[(long) (pow(2, j) - 1)], n)) {
                    rank++;
                    if (rank == n) {
                        t = cpuTime() - t1;
                        cout << "finished!" << endl;
                        glp_delete_prob(lp);
                        return;
                    }
                    rowechform(Arref, J, A[(long) (pow(2, j) - 1)], n, rank);
                    Asettled[rank - 1] = A[(long) (pow(2, j) - 1)];
                    settled_values[rank - 1] = v[(long) (pow(2, j) - 1)];
                }
                unsettled[(long) (pow(2, j) - 1)] = false;
                unsettled[(long) (s - pow(2, j))] = false;
                unsettled_p[j] = false;
            }
        }
    }

    glp_delete_prob(lp);

    double xS;
    while (rank < n) {
        iteration(unsettled, s, xS, n, A, x, v, epsi, prec, Arref, J, rank, Asettled, settled_values, iter, sr,
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
               vector<bool> &J, int &rank, vector<vector<bool>> &Asettled, vector<double> &settled_values,
               unsigned short int &iter, unsigned int &sr, vector<bool> &unsettled_p, vector<double> &singleton_bounds,
               bool &nlsu) {

    int t_size = 0;
    int t2_size = 0;
    vector<bool> T(s, false);
    vector<bool> T2(n, false);
    vector<unsigned int> T_coord(0, 0);
    vector<unsigned int> T2_coord(0, 0);

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
                T2_coord.push_back((long) (pow(2, i) - 1));
                t2_size++;
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
    subroutine(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, Asettled, sr, settled_values,
               unsettled, T_coord, s, epsi, v, T2_coord);
    if (rank == n) {
        return;
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
        if (unsettled_p[i] && !unsettled[(long) (pow(2, i) - 1)]) {
            unsettled_p[i] = false;
        }
    }

    // GLPK
    glp_prob *lp;
    int ia[1 + 1000] = {0};
    int ja[1 + 1000] = {0};
    double ar[1 + 1000] = {0};
    lp = glp_create_prob();
    glp_set_prob_name(lp, "P(2)");
    glp_set_obj_dir(lp, GLP_MIN);

    // objective function
    glp_add_cols(lp, n + 1);

    int ia_index = 1;

    for (unsigned short int j = 0; j < n; j++) {
        // set column names and coefficients
        const int col_index = j + 1;
        const std::string col_name_string = "x{" + std::to_string(col_index) + "}";
        const char *col_name = col_name_string.c_str();
        glp_set_col_name(lp, col_index, col_name);
        glp_set_col_bnds(lp, col_index, GLP_FR, 0.0, 0.0);
        glp_set_obj_coef(lp, col_index, 0.0);
    }

    // the only coefficient of the objective function, epsi
    int epsi_index = n + 1;
    glp_set_col_name(lp, epsi_index, "e(1)");
    glp_set_col_bnds(lp, epsi_index, GLP_FR, 0.0, 0.0);
    glp_set_obj_coef(lp, epsi_index, 1);

    vector<int> unsett_ineq_indices(s, -1);
    vector<int> impu_constr_indices(n, -1);
    int constr_index = 1;

    for (unsigned short int j = 0; j < n; j++) {
        if (unsettled_p[j]) {
            const int col_index = j + 1;

            // x({i}) >= v({i}); for all i = 1 ... n;
            constr_index++;
            glp_add_rows(lp, 1);
            impu_constr_indices[j] = constr_index;
            const std::string row_name_string =
                    "x({" + std::to_string(col_index) + "}) >= v({" + std::to_string(j) + "})" + " = " +
                    std::to_string(v[j]);
            const char *row_name = row_name_string.c_str();
            glp_set_row_name(lp, constr_index, row_name);
            glp_set_row_bnds(lp, constr_index, GLP_LO, singleton_bounds[j], 0.0);
            ia[ia_index] = constr_index;
            ja[ia_index] = col_index;
            ar[ia_index] = 1.0;
            ia_index++;
        }
    }

    // x(N) = v(N)
    constr_index++;
    glp_add_rows(lp, 1);
    const char *row_name = "x(N) = v(N)";
    glp_set_row_name(lp, constr_index, row_name);
    glp_set_row_bnds(lp, constr_index, GLP_FX, v[s], 0.0);
    for (unsigned short int j = 0; j < n; j++) {
        const int col_index = j + 1;
        ia[ia_index] = constr_index;
        ja[ia_index] = col_index;
        ar[ia_index] = 1.0;
        ia_index++;
    }

    for (unsigned int i = 0; i < s; i++) {
        if (unsettled[i]) {
            // x(S) + e(k) >= v(S); S in unsettled;
            constr_index++;
            glp_add_rows(lp, 1);
            unsett_ineq_indices[i] = constr_index;
            glp_set_row_name(lp, constr_index, "x(S) + e(k) >= v(S)");
            glp_set_row_bnds(lp, constr_index, GLP_LO, v[i], 0.0);
            ia[ia_index] = constr_index;
            ja[ia_index] = epsi_index;
            ar[ia_index] = 1.0;
            ia_index++;

            for (unsigned short int j = 0; j < n; j++) {
                const int col_index = j + 1;
                if (A[i][j]) {
                    ia[ia_index] = constr_index;
                    ja[ia_index] = col_index;
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
        glp_set_row_name(lp, constr_index, "x(S) + e(1) >= v(S)");
        glp_set_row_bnds(lp, constr_index, GLP_FX, settled_values[i], 0.0);

        for (unsigned short int j = 0; j < n; j++) {
            const int col_index = j + 1;
            if (Asettled[i][j]) {
                ia[ia_index] = constr_index;
                ja[ia_index] = col_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
        }
    }

    // NO warm start option. Kept going with solving from scratch

    glp_load_matrix(lp, ia_index - 1, ia, ja, ar);
    glp_simplex(lp, nullptr);
    iter++;

    for (unsigned short int j = 0; j < n; j++) {
        const int col_index = j + 1;
        x[j] = glp_get_col_prim(lp, col_index);
    }

    epsi = glp_get_col_prim(lp, epsi_index);

    for (unsigned int i = 0; i < s; i++) {
        if (unsettled[i]) {
            if (glp_get_row_dual(lp, unsett_ineq_indices[i]) > prec) {
                if (binrank(Arref, J, A[i], n)) {
                    rank++;
                    if (rank == n) {
                        glp_delete_prob(lp);
                        return;
                    }
                    rowechform(Arref, J, A[i], n, rank);
                    Asettled[rank - 1] = A[i];
                    settled_values[rank - 1] = v[i] - epsi;
                }
                unsettled[i] = false;
                unsettled[s - 1 - i] = false;
            }
        }
    }

    for (unsigned short int j = 0; j < n; j++) {
        if (unsettled_p[j]) {
            if (glp_get_row_dual(lp, impu_constr_indices[j]) > prec) {
                if (binrank(Arref, J, A[(long) (pow(2, j) - 1)], n)) {
                    rank++;
                    if (rank == n) {
                        glp_delete_prob(lp);
                        return;
                    }
                    rowechform(Arref, J, A[(long) (pow(2, j) - 1)], n, rank);
                    Asettled[rank - 1] = A[(long) (pow(2, j) - 1)];
                    settled_values[rank - 1] = v[(long) (pow(2, j) - 1)];
                }
                unsettled[(long) (pow(2, j) - 1)] = false;
                unsettled[(long) (s - pow(2, j))] = false;
                unsettled_p[j] = false;
            }
        }
    }

    glp_delete_prob(lp);
}

void subroutine(vector<bool> &U, vector<bool> &U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2,
                vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n, int &t_size,
                int &t2_size, int &rank, vector<vector<bool>> &Asettled, unsigned int &sr,
                vector<double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s,
                double &epsi, vector<double> &v, vector<unsigned int> &T2_coord) {
    // GLPK
    glp_prob *lp;
    int ia[1 + 1000] = {0};
    int ja[1 + 1000] = {0};
    double ar[1 + 1000] = {0};
    lp = glp_create_prob();
    glp_set_prob_name(lp, "SUBROUTINE(k)");
    glp_set_obj_dir(lp, GLP_MAX);

    // objective function
    glp_add_cols(lp, t_size + t2_size + rank);

    for (int i = 0; i < t_size + t2_size; i++) {
        // sr_obj from 0 -> t_size + t2_size
        const int col_index = i + 1;
        const std::string col_name_string = "sr{" + std::to_string(col_index) + "}";
        const char *col_name = col_name_string.c_str();
        glp_set_col_name(lp, col_index, col_name);
        glp_set_col_bnds(lp, col_index, GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(lp, col_index, 1);
    }

    for (int i = t_size + t2_size; i < t_size + t2_size + rank; i++) {
        // sr_obj from t_size + t2_size -> rank + t_size + t2_size
        const int col_index = i + 1;
        const std::string col_name_string = "sr{" + std::to_string(col_index) + "}";
        const char *col_name = col_name_string.c_str();
        glp_set_col_name(lp, col_index, col_name);
        glp_set_col_bnds(lp, col_index, GLP_FR, 0.0, 0.0);
        glp_set_obj_coef(lp, col_index, 0);
    }

    vector<vector<int>> bal_indices(t_size + t2_size + rank, vector<int>(n + 1, -1));
    vector<int> constr_indices(n + 1, -1);
    int ia_index = 1;
    int constr_index = 0;

    for (unsigned int j = 0; j < n; j++) {
        // eq == 0
        constr_index++;
        constr_indices[j] = constr_index;
        glp_add_rows(lp, 1);
        const std::string row_name_string = "t({" + std::to_string(j + 1) + "}) = 0";
        const char *row_name = row_name_string.c_str();
        glp_set_row_name(lp, constr_index, row_name);
        glp_set_row_bnds(lp, constr_index, GLP_FX, 0, 0.0);

        for (int i = 0; i < t_size; i++) {
            if (Atight[i][j]) {
                int col_index = i + 1;
                bal_indices[i][j] = ia_index;
                ia[ia_index] = constr_index;
                ja[ia_index] = col_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
            if (rank > i && Asettled[i][j]) {
                int col_index = i + t_size + t2_size + 1;
                bal_indices[i][j] = ia_index;
                ia[ia_index] = constr_index;
                ja[ia_index] = col_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
        }
        for (int i = 0; i < t2_size; i++) {
            if (Atight2[i][j]) {
                int col_index = i + t_size + 1;
                bal_indices[i][j] = ia_index;
                ia[ia_index] = constr_index;
                ja[ia_index] = col_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
        }
        for (int i = t_size; i < rank; i++) {
            if (rank > t_size && Asettled[i][j]) {
                int col_index = i + t_size + t2_size + 1;
                bal_indices[i][j] = ia_index;
                ia[ia_index] = constr_index;
                ja[ia_index] = col_index;
                ar[ia_index] = 1.0;
                ia_index++;
            }
        }
    }

    // pos_eq == 1
    constr_index++;
    constr_indices[n] = constr_index;
    glp_add_rows(lp, 1);
    const std::string row_name_string = "t({" + std::to_string(n) + "}) = 1";
    const char *row_name = row_name_string.c_str();
    glp_set_row_name(lp, constr_index, row_name);
    glp_set_row_bnds(lp, constr_index, GLP_FX, 1, 0.0);

    for (int i = 0; i < t_size; i++) {
        const int col_index = i + 1;
        bal_indices[i][n] = ia_index;
        ia[ia_index] = constr_index;
        ja[ia_index] = col_index;
        ar[ia_index] = 1.0;
        ia_index++;
    }

    glp_load_matrix(lp, ia_index - 1, ia, ja, ar);
    int feas = glp_simplex(lp, nullptr);

    sr++;

    vector<double> u(t_size + t2_size, 0);

    if (feas == 0) {
        for (int i = 0; i < t_size + t2_size; i++) {
            const int col_index = i + 1;
            u[i] = glp_get_col_prim(lp, col_index);
        }
    }

    unsigned int sumt = 0;
    unsigned short int sumt2 = 0;
    vector<bool> t(t_size, false);
    vector<bool> t2(t2_size, false);

    while (feas) {
        subr_upd(Arref, J, lp, ia_index, constr_indices, bal_indices, ia, ja, ar, n, prec, U, U2, sumt, sumt2, t, t2,
                 Atight, Atight2, t_size, t2_size, rank, Asettled, settled_values, unsettled, T_coord, s, epsi, v,
                 T2_coord, u);
        if (rank == n) {
            return;
        } else {
            int i = 0;
            while (i < t_size) {
                if (!t[i] && !(binrank(Arref, J, Atight[i], n))) {
                    U[i] = false;
                    t[i] = true;

                    int ia_ref_index = bal_indices[i][n];
                    if (ia_ref_index > -1) {
                        ar[ia_ref_index] = ar[ia_ref_index] - 1;
                    } else {
                        const int constr_ref_index = constr_indices[n];
                        const int col_index = i + 1;
                        bal_indices[i][n] = ia_index;
                        ia[ia_index] = constr_ref_index;
                        ja[ia_index] = col_index;
                        ar[ia_index] = -1;
                        ia_index++;
                    }

                    const int sr_obj_index = i + 1;
                    const double sr_obj_coef = glp_get_obj_coef(lp, sr_obj_index);
                    glp_set_obj_coef(lp, sr_obj_index, sr_obj_coef - 1);

                    for (unsigned int j = 0; j < n; j++) {
                        if (Atight[i][j]) {
                            ia_ref_index = bal_indices[i][j];
                            if (ia_ref_index > -1) {
                                ar[ia_ref_index] = ar[ia_ref_index] - 1;
                            } else {
                                const int constr_ref_index = constr_indices[j];
                                const int col_index = i + 1;
                                bal_indices[i][n] = ia_index;
                                ia[ia_index] = constr_ref_index;
                                ja[ia_index] = col_index;
                                ar[ia_index] = -1;
                                ia_index++;
                            }
                        }
                    }
                    sumt++;
                    unsettled[T_coord[i]] = false;
                    unsettled[s - 1 - T_coord[i]] = false;
                }
                i++;
            }
            i = 0;
            while (i < t2_size) {
                if (!t2[i]) {
                    if (!(binrank(Arref, J, Atight2[i], n))) {
                        U2[i] = false;
                        t2[i] = true;

                        const int sr_obj_index = i + t_size + 1;
                        const double sr_obj_coef = glp_get_obj_coef(lp, sr_obj_index);
                        glp_set_obj_coef(lp, sr_obj_index, sr_obj_coef - 1);

                        for (unsigned short int j = 0; j < n; j++) {
                            if (Atight2[i][j]) {
                                const int ia_ref_index = bal_indices[i + t_size][j];
                                if (ia_ref_index > -1) {
                                    ar[ia_ref_index] = ar[ia_ref_index] - 1;
                                } else {
                                    const int constr_ref_index = constr_indices[j];
                                    const int col_index = i + t_size + 1;
                                    bal_indices[i][n] = ia_index;
                                    ia[ia_index] = constr_ref_index;
                                    ja[ia_index] = col_index;
                                    ar[ia_index] = -1;
                                    ia_index++;
                                }
                            }
                        }
                        sumt2++;
                        if (unsettled[T2_coord[i]]) {
                            unsettled[T2_coord[i]] = false;
                        }
                        if (unsettled[s - 1 - T2_coord[i]]) {
                            unsettled[s - 1 - T2_coord[i]] = false;
                        }
                    }
                }
                i++;
            }
        }
        if (sumt == t_size) {
            return;
        } else {
            glp_load_matrix(lp, ia_index - 1, ia, ja, ar);
            feas = glp_simplex(lp, nullptr);

            sr++;
            if (feas == 0) {
                for (int i = 0; i < t_size + t2_size; i++) {
                    const int col_index = i + 1;
                    u[i] = glp_get_col_prim(lp, col_index);
                }
            }
        }
    }
    glp_delete_prob(lp);
}

void subr_upd(vector<vector<double>> &Arref, vector<bool> &J, glp_prob *lp, int &ia_index, vector<int> constr_indices,
              vector<vector<int>> bal_indices, int ia[], int ja[], double ar[], unsigned short int &n, double &prec,
              vector<bool> &U, vector<bool> &U2, unsigned int &sumt, unsigned short int &sumt2, vector<bool> &t,
              vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, int &t_size, int &t2_size,
              int &rank, vector<vector<bool>> &Asettled, vector<double> &settled_values, vector<bool> &unsettled,
              vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector<double> &v,
              vector<unsigned int> &T2_coord, vector<double> &u) {

    int i = 0;
    while (i < t_size && sumt < t_size) {
        if (!t[i] && u[i] > prec) {
            U[i] = false;
            t[i] = true;

            int ia_ref_index = bal_indices[i][n];
            if (ia_ref_index > -1) {
                ar[ia_ref_index] = ar[ia_ref_index] - 1;
            } else {
                const int constr_ref_index = constr_indices[n];
                const int col_index = i + 1;
                bal_indices[i][n] = ia_index;
                ia[ia_index] = constr_ref_index;
                ja[ia_index] = col_index;
                ar[ia_index] = -1;
                ia_index++;
            }

            sumt++;
            if (binrank(Arref, J, Atight[i], n)) {
                if (rank == n - 1) {
                    rank++;
                    return;
                }
                rowechform(Arref, J, Atight[i], n, rank);
                rank++;
                Asettled[rank - 1] = Atight[i];
                settled_values[rank - 1] = v[T_coord[i]] - epsi;
                unsettled[T_coord[i]] = false;
                unsettled[s - 1 - T_coord[i]] = false;

                const int col_index = i + 1;
                glp_set_col_bnds(lp, col_index, GLP_FR, 0.0, 0.0);
            } else {
                for (unsigned short int j = 0; j < n; j++) {
                    if (Atight[i][j]) {
                        ia_ref_index = bal_indices[i][j];
                        if (ia_ref_index > -1) {
                            ar[ia_ref_index] = ar[ia_ref_index] - 1;
                        } else {
                            const int constr_ref_index = constr_indices[j];
                            const int col_index = i + 1;
                            bal_indices[i][n] = ia_index;
                            ia[ia_index] = constr_ref_index;
                            ja[ia_index] = col_index;
                            ar[ia_index] = -1;
                            ia_index++;
                        }
                    }
                }
                unsettled[T_coord[i]] = false;
                unsettled[s - 1 - T_coord[i]] = false;
            }
        }
        i++;
    }
    i = 0;
    while (i < t2_size && sumt2 < t2_size) {
        if (!t2[i] && u[i + t_size] > prec) {
            U2[i] = false;
            t2[i] = true;
            sumt2++;

            const int sr_obj_index = i + t_size + 1;
            const double sr_obj_coef = glp_get_obj_coef(lp, sr_obj_index);
            glp_set_obj_coef(lp, sr_obj_index, sr_obj_coef - 1);

            if (binrank(Arref, J, Atight2[i], n)) {
                if (rank == n - 1) {
                    rank++;
                    return;
                }
                rowechform(Arref, J, Atight2[i], n, rank);
                rank++;
                Asettled[rank - 1] = Atight2[i];
                settled_values[rank - 1] = v[T2_coord[i]];
                unsettled[T2_coord[i]] = false;
                unsettled[s - 1 - T2_coord[i]] = false;

                const int col_index = i + t_size + 1;
                glp_set_col_bnds(lp, col_index, GLP_FR, 0.0, 0.0);
            } else {
                for (unsigned short int j = 0; j < n; j++) {
                    if (Atight2[i][j]) {
                        const int ia_ref_index = bal_indices[i + t_size][j];
                        if (ia_ref_index > -1) {
                            ar[ia_ref_index] = ar[ia_ref_index] - 1;
                        } else {
                            const int constr_ref_index = constr_indices[j];
                            const int col_index = i + t_size + 1;
                            bal_indices[i][n] = ia_index;
                            ia[ia_index] = constr_ref_index;
                            ja[ia_index] = col_index;
                            ar[ia_index] = -1;
                            ia_index++;
                        }
                    }
                }
                unsettled[T2_coord[i]] = false;
                unsettled[s - 1 - T2_coord[i]] = false;
            }
        }
        i++;
    }
}