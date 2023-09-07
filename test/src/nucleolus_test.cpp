//
// Created by secci on 25/07/2023.
//
#include "cgt_test.h"

#define TEST_FILES_DIR _TEST_FILES_DIR

TEST_CASE("Nucleolus")
{
    std::string tests_file = TEST_FILES_DIR;
    tests_file.append("/games.json");
    ifstream ifs(tests_file);
    Json::Reader reader;
    Json::Value games;
    reader.parse(ifs, games);

    for (int i = 0; i < games.size(); i++)
    {
        unsigned short int n = games[i]["n"].asUInt();

        const std::string &test_name = games[i]["name"].asString();

        std::vector<double> v;
        const Json::Value &v_in = games[i]["v"];
        for (int j = 0; j < v_in.size(); j++)
        {
            v.push_back(v_in[j].asDouble());
        }

        std::vector<double> ground_nucleolus;
        const Json::Value &nucleolus_in = games[i]["nucleolus"];
        for (int j = 0; j < nucleolus_in.size(); j++)
        {
            ground_nucleolus.push_back(nucleolus_in[j].asDouble());
        }

        SUBCASE(test_name.c_str())
        {
            unsigned int s = pow(2, n) - 2;
            std::vector<double> x(n, 0);
            std::vector<double> singleton_bounds(n, 0);
            std::vector<double> excess(s, 0);
            std::vector<bool> unsettled(s + 1, true);
            unsettled[s] = false;

            bool disp = false;
            bool nlsu = false;
            double impu = 0;
            double prec = pow(10, -3);
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
            std::vector<std::vector<bool>> A(s + 1, std::vector<bool>(n, false));
            A_mx(A, n, s);
            excess_init(excess, unsettled, A, x, v, s, n);
            nucleolus(disp, n, s, excess, prec, unsettled, iter, piv, sr, t, x, A, t1, singleton_bounds, nlsu);

            for (unsigned int i = 0; i < n; i++)
            {
                CHECK(abs(x[i] - ground_nucleolus[i]) <= prec);
            }
        }
    }
}
