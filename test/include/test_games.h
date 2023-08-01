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
     }, {
        3,
        "n = 3 game on electricity by SatyaRamesh and Radhakrishna (2009)",
        {1.275, 3.471, 7.005, 1.466, 4.081, 5.672, 11.210},
        {3.0893, 5.2853, 2.8355}
     }};
#endif