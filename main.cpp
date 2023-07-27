#include "lib/include/primal_dual.h"

int main() {
    unsigned short int n = 3;
    double x[] = {0, 0, 0};
    double v[] = {0, 0, 0, 5, 5, 1, 9};
    primal_dual(v, x, n);
    for (unsigned int i = 0; i < n; i++)
        cout << fixed << setprecision(17) << x[i] << endl;
    return 0;
}