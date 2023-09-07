//
// Created by secci on 25/07/2023.
//
#include "cgt_test.h"

#define TEST_FILES_DIR _TEST_FILES_DIR

TEST_CASE("Shapley")
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

        std::vector<double> ground_shapley;
        const Json::Value &shapley_in = games[i]["shapley"];
        for (int j = 0; j < shapley_in.size(); j++)
        {
            ground_shapley.push_back(shapley_in[j].asDouble());
        }

        SUBCASE(test_name.c_str())
        {
            double prec = pow(10, -3);

            vector<const shapley::OrdinalPlayer *> players;
            for (unsigned short int i = 0; i < n; i++)
            {
                players.push_back(new shapley::OrdinalPlayer(i, v));
            }

            shapley::OrdinalCharacteristicFunction char_func(v);

            map<unsigned int short, double> shapley_values_map = compute(players, char_func);
            vector<double> shapley_values_vec;

            for (auto elem : shapley_values_map)
            {
                shapley_values_vec.push_back(elem.second);
            }

            for (unsigned int i = 0; i < n; i++)
            {
                CHECK(abs(shapley_values_vec[i] - ground_shapley[i]) <= prec);
            }

            for (unsigned short int i = 0; i < n; i++)
            {
                delete players[i];
            }
        }
    }
}
