#include "bnf_module.h"

using namespace emscripten;

vector<double> bnf_run(vector<double> v_in, unsigned short int n_in)
{
    unsigned short int n = n_in;
    unsigned int s = pow(2, n) - 2;
    vector<double> x(n, 0);
    vector<double> singleton_bounds(n, 0);
    vector<double> v(v_in);
    vector<double> excess(s, 0);
    vector<bool> unsettled(s + 1, true);
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
    vector<vector<bool>> A(s + 1, vector<bool>(n, false));
    A_mx(A, n, s);
    excess_init(excess, unsettled, A, x, v, s, n);
    bnf(disp, n, s, excess, prec, unsettled, iter, piv, sr, t, x, A, t1, singleton_bounds, nlsu);
    return x;
}

EMSCRIPTEN_BINDINGS(bnf_module)
{
    register_vector<double>("DoubleList");
    emscripten::function("bnf_run", &bnf_run);
}