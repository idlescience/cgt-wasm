//
// Created by secci on 25/07/2023.
//
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "nucleolus_test.h"

TEST_CASE("Primal-Dual")
{
    SUBCASE("should work in a simple 3 players game")
    {
        unsigned short int n = 3;
        double x[] = {0, 0, 0};
        double v[] = {0, 0, 0, 3, 0, 1, 4};
        primal_dual(v, x, n);
        double x_0 = round(x[0] * 100.0) / 100.0;
        CHECK(x_0 == 2.50);
    }
}

TEST_CASE("BNF")
{
    for (game test_game : GAMES)
    {
        const std::string test_name = "should work in a " + test_game.name + " game";
        SUBCASE(test_name.c_str())
        {
            srand(test_game.seed);
            unsigned short int n = test_game.n;
            unsigned short int type = test_game.type;
            unsigned int s = pow(2, n) - 2;
            vector<double> x(n, 0);
            vector<double> singleton_bounds(n, 0);
            vector<double> v(s + 1, 0);
            vector<double> excess(s, 0);
            vector<bool> unsettled(s + 1, true);
            unsettled[s] = false;
            cout << "Generating game...";
            if (type == 1)
            {
                type1(v, s, n);
            }
            else if (type == 2)
            {
                type2(v, s, n);
            }
            else if (type == 4)
            {
                type4(v, s, n);
            }
            cout << "done!" << endl;
            cout << "Running BNF..." << endl;

            bool disp = true;
            bool nlsu = false;
            double impu = 0;
            double prec = pow(10, -6);
            unsigned short int iter = 0;
            unsigned int piv = 0;
            unsigned int sr = 0;
            double t = 0;
            double t1 = cpuTime();

            for (unsigned short int i = 0; i < n; i++)
            {
                singleton_bounds[i] = v[pow(2, i) - 1];
                impu += singleton_bounds[i];
            }
            x = singleton_bounds;
            for (unsigned short int i = 0; i < n; i++)
            {
                x[i] += (v[s] - impu) / n;
            }
            vector<vector<bool>> A(s + 1, vector<bool>(n, false));
            A_mx(A, n, s);
            excess_init(excess, unsettled, A, x, v, s, n);
            // bnf(disp, n, s, excess, prec, unsettled, iter, piv, sr, t, x, A, t1, singleton_bounds, nlsu);

            for (unsigned int i = 0; i < n; i++)
            {
                double x_rounded = round(x[i] * 100.0) / 100.0;
                double payoff_rounded = round(test_game.payoff[i] * 100.0) / 100.0;
                CHECK(x_rounded == payoff_rounded);
            }
        }
    }
}
