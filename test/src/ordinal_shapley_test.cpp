//
// Created by secci on 25/07/2023.
//
#include "cgt_test.h"
using namespace Shapley;

TEST_CASE("Ordinal Shapley")
{
    for (auto test_game : GAMES)
    {
        const std::string test_name = "should work on " + test_game.name;
        SUBCASE(test_name.c_str())
        {
            unsigned short int n_in = test_game.n;
            vector<double> v_in(test_game.v);
            double prec = pow(10, -3);

            vector<const OrdinalPlayer *> players;
            for (unsigned short int i = 0; i < n_in; i++)
            {
                players.push_back(new OrdinalPlayer(i, v_in));
            }

            OrdinalCharacteristicFunction char_func(v_in);

            map<const OrdinalPlayer *, double> shapley_values_map = compute(players, char_func);
            vector<double> shapley_values_vec;

            for (auto elem : shapley_values_map)
            {
                shapley_values_vec.push_back(elem.second);
            }

            for (unsigned short int i = 0; i < n_in; i++)
            {
                delete players[i];
            }

            for (unsigned int i = 0; i < n_in; i++)
            {
                CHECK(abs(shapley_values_vec[i] - test_game.shapley[i]) <= prec);
            }
        }
    }
}
