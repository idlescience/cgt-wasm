#include "lib/nucleolus.h"

int main() {
    unsigned short int n = 3;
    double x[] = {0, 0, 0};
    double v[] = {0, 0, 0, 5, 5, 1, 9};
    PD(v, x, n);
    for (unsigned int i = 0; i < n; i++)
        cout << fixed << setprecision(17) << x[i] << endl;
    cout << "Press 0 then Enter to quit: ";
    double quit;
    cin >> quit;
    cin.get();
    return 0;
}