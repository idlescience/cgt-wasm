#ifndef _TEST_GAMES_H
#define _TEST_GAMES_H
struct game {
    int n;
    std::string name;
    vector<double> v;
    vector<double> nucleolus;
};
static const game GAMES[] = {{
        3,
        "n = 3 game Reijinierse (1995)",
        {0, 0, 5, 0, 5, 1, 9},
        {5, 2, 2}
    //  }, {
    //     3,
    //     "n = 3 game on electricity by SatyaRamesh and Radhakrishna (2009)",
    //     {1.275, 3.471, 1.466, 7.005, 4.081, 5.672, 11.210},
    //     {3.0893, 5.2853, 2.8355}
    //  }, {
    //     3,
    //     "n = 3 game on manufacturing by Lemaire (1991)",
    //     {46125.0, 17437.5, 5812.5, 69187.5, 53812.5, 30750.0, 90000.0},
    //     {52687.5, 24468.8, 12843.8}
    // }, {
    //     3,
    //     "n = 3 game on manufacturing by Oh and Shin (2012)",
    //     {-375144, -245280, -211239, -568232, -551055, -452411, -772500},
    //     {-331695.3, -233051.3, -207753.5}
    }};
#endif