//
// Created by secci on 25/07/2023.
//
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "bnf.h"
#include "primal_dual.h"
#include "nucleolus_test.h"

TEST_CASE("Simple primal_dual of 3 players") {
    unsigned short int n = 3;
    double x[] = {0, 0, 0};
    double v[] = {0, 0, 0, 5, 5, 1, 9};
    primal_dual(v, x, n);
    double x_0 = round(x[0] * 100.0) / 100.0;
    CHECK(x_0 == 1.33);
}


TEST_CASE("RStudio primal_dual of 3 players") {
    unsigned short int n = 3;
    double x[] = {0, 0, 0};
    double v[] = {0, 0, 0, 3, 0, 1, 4};
    primal_dual(v, x, n);
    double x_0 = round(x[0] * 100.0) / 100.0;
    CHECK(x_0 == 2.50);
}


TEST_CASE("Benedek BNF test for type 1 game") {
    unsigned short int n = 0;
    unsigned short int type = 0;
    unsigned int seed = 0;
    bool disp = false;
    bool memo = false;
    bool nlsu = false;
    ifstream inp;
    cout << "Reading the input...";
    inp.open("input.txt");
    inp >> n >> type >> seed >> disp >> memo >> nlsu;
    inp.close();
    cout << "done!" << endl;
    if (seed == 0) {
        seed = GetTickCount();
    }
    srand(seed);
    unsigned int s = pow(2, n) - 2;
    vector<double> x(n, 0);
    vector<double> singleton_bounds(n, 0);
    double impu = 0;
    double prec = pow(10, -6);
    vector<double> excess(s, 0);
    vector<bool> unsettled(s + 1, true);
    unsettled[s] = false;
    unsigned short int iter = 0;
    unsigned int piv = 0;
    unsigned int sr = 0;
    double t = 0;

    vector<double> v(s + 1, 0);

    double t1 = cpuTime();
    for (unsigned short int i = 0; i < n; i++) {
        singleton_bounds[i] = v[pow(2, i) - 1];
        impu += singleton_bounds[i];
    }
    x = singleton_bounds;
    for (unsigned short int i = 0; i < n; i++) {
        x[i] += (v[s] - impu) / n;
    }

    vector<vector<bool>> A(s + 1, vector<bool>(n, false));
    A_mx(A, n, s);
    excess_init(excess, unsettled, A, x, v, s, n);
    BNF(disp, n, s, excess, prec, unsettled, iter, piv, sr, t, x, A, t1, singleton_bounds, nlsu);

    ofstream res;
    res.open("results.txt", ofstream::out | ofstream::trunc);
    res << seed << endl << t << endl << iter << endl << piv << endl << sr << endl;
    for (unsigned int i = 0; i < n; i++) {
        res << fixed << setprecision(17) << x[i] << endl;
    }
    res.close();
    cout << "Press 0 then Enter to quit: ";
    double quit;
    cin >> quit;
    cin.get();
    double x_0 = round(x[0] * 100.0) / 100.0;
    CHECK(x_0 == 2.50);
}
