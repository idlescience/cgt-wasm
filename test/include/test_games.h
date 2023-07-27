#ifndef _TEST_GAMES_H
#define _TEST_GAMES_H
struct game {
    int seed;
    int n;
    int type;
    vector<double> payoff;
};
static const game GAMES[] = {{461125, 5, 1,
                              {60.00000000000002132, 112.99999999999998579, 37.50000000000000000, 145.50000000000000000,
                               106.99999999999998579}}};
#endif