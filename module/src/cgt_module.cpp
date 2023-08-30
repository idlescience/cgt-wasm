#include "cgt_module.h"

using namespace emscripten;

vector<double> nucleolus_run(vector<double> v_in, unsigned short int n_in)
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
    nucleolus(disp, n, s, excess, prec, unsettled, iter, piv, sr, t, x, A, t1, singleton_bounds, nlsu);
    return x;
}

vector<double> shapley_run(vector<double> v_in, unsigned short int n_in)
{
    using namespace Shapley;

    class ModulePlayer : public Player
    {
    public:
        explicit ModulePlayer(int index, const vector<double>& v) : position(index), v_ref(v) {}
        double getContribution() const override
        {
            return v_ref[1 << position];
        }

    protected:
        int position;
        const vector<double>& v_ref;
    };

    class ModuleCharacteristicFunction : public CharacteristicFunction<ModulePlayer>
    {
    public:
        double getValue(const Coalition<ModulePlayer> &coalition) const override
        {
            double value = 0.0;
            // Find largest taxi fare of any coalition member.
            for (const ModulePlayer *member : coalition.getMembers())
            {
                if (member->getContribution() > value)
                    value = member->getContribution();
            }
            return value;
        }
    };
}

EMSCRIPTEN_BINDINGS(cgt)
{
    register_vector<double>("DoubleVector");
    emscripten::function("nucleolus", &nucleolus_run);
    emscripten::function("shapley", &shapley_run);
}