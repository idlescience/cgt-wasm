//
// Created by secci on 25/07/2023.
//
#include "cgt_test.h"
using namespace Shapley;

TEST_CASE("Power Set")
{
    const game test_game = GAMES[0];
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
        Coalition<OrdinalPlayer> grand_coalition(players);

        std::vector<Coalition<OrdinalPlayer>> ans = powerSet(grand_coalition);

        CHECK(ans.size() == pow(2, n_in) - 1);
        CHECK(ans.at(0).getPlayerAt(0)->getPosition() == 0);
        CHECK(ans.at(10).getPlayerAt(0)->getPosition() == 0);
        CHECK(ans.at(10).getPlayerAt(1)->getPosition() == 1);
    }
}

TEST_CASE("Power Set")
{
    const game test_game = GAMES[0];
    const std::string test_name = "should display lexicographically ordened payoffs" + test_game.name;
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
        Coalition<OrdinalPlayer> grand_coalition(players);
        OrdinalCharacteristicFunction char_func(v_in);

        std::vector<Coalition<OrdinalPlayer>> ans = powerSet(grand_coalition);

        ofstream myfile;
        myfile.open("r_coalition.txt");
        for (unsigned int i = 0; i < ans.size(); i++)
        {
            Coalition<OrdinalPlayer> coalition = ans.at(i);
            const double value = char_func.getValue(coalition);
            myfile << value << "\n";
        }
        myfile.close();
    }
}

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

            for (unsigned int i = 0; i < n_in; i++)
            {
                CHECK(abs(shapley_values_vec[i] - test_game.shapley[i]) <= prec);
            }

            for (unsigned short int i = 0; i < n_in; i++)
            {
                delete players[i];
            }
        }
    }
}
