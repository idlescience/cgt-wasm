#include "PD.h"

int main() {
    unsigned short int n = 3;
    unsigned int seed = 0;
    bool disp = true;
    bool nlsu = false;
    seed = GetTickCount();
    srand(seed);
    unsigned int s = pow(2, n) - 2;
    vector<double> x(n, 0);
    vector<double> v{0, 0, 0, 5, 5, 1, 9};
    unsigned short int iter = 0;
    unsigned int piv = 0;
    unsigned int sr = 0;
    double t = 0;
    PD(disp, n, v, iter, piv, sr, t, x, s, nlsu);
    cout << seed << endl << t << endl << iter << endl << piv << endl;
    for (unsigned int i = 0; i < n; i++)
        cout << fixed << setprecision(17) << x[i] << endl;
    cout << "Press 0 then Enter to quit: ";
    double quit;
    cin >> quit;
    cin.get();
    return 0;
}