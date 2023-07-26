#include "common.h"

double cpuTime() {
    return (double) clock() / CLOCKS_PER_SEC;
}

void rowechform(vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &B, unsigned short int &n, int &rank) {
    double prec = pow(10, -10);
    vector<vector<double>> rref(rank + 1, vector<double>(n, 0));
    rref[0] = Arref[0];// first row done
    for (unsigned int i = 1; i < rank + 1; i++) {
        for (unsigned short int j = 0; j < n; j++) {
            if (i < rank) {
                if (Arref[i][j] > prec || Arref[i][j] < -prec)
                    rref[i][j] = Arref[i][j];
            } else {
                if (B[j])
                    rref[i][j] = 1;
            }
        }
    }
    for (unsigned int i = 1; i < rank + 1; i++) {
        if (rref[i][0] > prec || rref[i][0] < -prec) {
            if (rref[i][0] < 1 + prec && rref[i][0] > 1 - prec)
                vec_subtract(rref[i], rref[0], rref[i]);
            else {
                rowechform_piv2(rref, i, n);
            }
        }
    }

    unsigned int k = 1;
    unsigned short int j = 1;
    while (k < rank + 1 && j < n) {
        rowechform_loop(rref, J, k, j, rank, prec, n);
    }
    if (rank + 1 < n) {
        for (int l = 1; l < rank + 1; l++) {
            Arref[l] = rref[l];
        }
    } else {
        for (int l = 1; l < n; l++) {
            Arref[l] = rref[l];
        }
    }
}

void rowechform_loop(vector<vector<double>> &rref, vector<bool> &J, unsigned int &i, unsigned short int &j, int &rank,
                     double &prec, unsigned short int &n) {
    vector<unsigned int> nonz(0, 0);
    vector<unsigned int> ones(0, 0);
    for (unsigned int k = i; k < rank + 1; k++) {
        if (rref[k][j] > prec || rref[k][j] < -prec) {
            if (rref[k][j] < 1 + prec && rref[k][j] > 1 - prec) {
                ones.push_back(k);
            } else {
                nonz.push_back(k);
            }
        }
    }
    if (ones.empty() && nonz.empty()) {
        j++;
    } else {
        if (ones.empty()) {
            if (nonz[0] != i) {
                swap_ith_and_firstnz(rref, nonz, i);// swap i-th and first non-zero row
            }
            sc_vec_prod(rref[i], 1 / rref[i][j], rref[i]);
        } else {
            if (ones[0] != i) {
                swap_ith_and_firstone(rref, ones, nonz, i);
            }
        }
        if (!ones.empty()) {
            if (ones[0] == i) {
                for (unsigned int k = 1; k < ones.size(); k++) {
                    vec_subtract(rref[ones[k]], rref[ones[k]], rref[i]);
                }
            } else {
                for (unsigned int one: ones) {
                    vec_subtract(rref[one], rref[one], rref[i]);
                }
            }
        }
        if (!nonz.empty()) {
            if (nonz[0] == i) {
                for (unsigned int k = 1; k < nonz.size(); k++) {
                    rowechform_piv(rref, nonz, i, j, k, n);
                }
            } else {
                for (unsigned int k = 0; k < nonz.size(); k++) {
                    rowechform_piv(rref, nonz, i, j, k, n);
                }
            }
        }
        i++;
        J[j] = false;
        j++;
    }
}

void rowechform_piv(vector<vector<double>> &rref, vector<unsigned int> &nonz, unsigned int &i, unsigned short int &j,
                    unsigned int &k, unsigned short int &n) {
    vector<double> aux(n, 0);
    sc_vec_prod(aux, rref[nonz[k]][j], rref[i]);
    vec_subtract(rref[nonz[k]], rref[nonz[k]], aux);
}

void rowechform_piv2(vector<vector<double>> &rref, unsigned int &i, unsigned short int &n) {
    vector<double> aux(n, 0);
    sc_vec_prod(aux, rref[i][0], rref[0]);
    vec_subtract(rref[i], aux, rref[i]);
}

void swap_ith_and_firstnz(vector<vector<double>> &rref, vector<unsigned int> &nonz, unsigned int &i) {
    vector<double> aux = rref[nonz[0]]; // swap i-th and first non-zero row
    rref[nonz[0]] = rref[i];
    rref[i] = aux;
    nonz[0] = i;
}

void swap_ith_and_firstone(vector<vector<double>> &rref, vector<unsigned int> &ones, vector<unsigned int> &nonz,
                           unsigned int &i) {
    vector<double> aux = rref[ones[0]];
    rref[ones[0]] = rref[i];
    rref[i] = aux;
    if (!nonz.empty()) {
        if (nonz[0] == i) {
            nonz[0] = ones[0];
        }
    }
    ones[0] = i;
}

bool binrank(vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &b, unsigned short int &n) {
    double prec = pow(10, -10);
    vector<double> B(n, 0);
    for (unsigned short int i = 0; i < n; i++) {
        if (b[i])
            B[i] = 1;
    }
    unsigned int m = 0;
    bool size = true;
    while (size) {
        if (nonz_vec(Arref[m], prec))
            m++;
        else
            size = false;
    }
    if (m >= n)
        return false;
    else {
        vector<bool> pivot_col(n, false);
        for (unsigned short int i = 0; i < n; i++) {
            if (!J[i]) {
                pivot_col[i] = true;
            }
        }
        unsigned short int j = 0;
        vector<bool> piv(n, false);
        vector<double> aux(n, 0);
        unsigned short int k = 0;
        unsigned short int I = 0;
        unsigned int s = 0;
        unsigned short int ind = 0;
        unsigned short int count = 0;
        while (j < n) {
            for (unsigned short i = 0; i < n; i++) {
                if (B[i] > prec || B[i] < -prec)
                    piv[i] = true;
            }
            sum_vecb(s, piv);
            if (s == 0)
                return false;
            else {
                while (k == 0) {
                    if (piv[I])
                        k = I + 1;
                    I++;
                }
                k--;
                I = 0;
                if (J[k])
                    return true;
                else {
                    while (count < k + 1) {
                        if (pivot_col[count])
                            ind++;
                        count++;
                    }
                    ind--;
                    count = 0;
                    sc_vec_prod(aux, B[k] / Arref[ind][k], Arref[ind]);
                    vec_subtract(B, B, aux);
                    j++;
                }
            }
            for (unsigned short int l = 0; l < n; l++)
                piv[l] = false;
            k = 0;
            ind = 0;
        }
        return false;
    }
}

void A_mx(vector<vector<bool>> &A, unsigned short int &n, unsigned int &s) {
    // creates boolean matrix A containing all the possible n-length boolean vectors (except for full zeros)
    for (unsigned int k = 0; k != s + 1; k++) {
        unsigned int i = 2;
        for (unsigned short int c = 0; c < n - 2; c++)
            i += i;
        unsigned int j = k + 1;
        unsigned short int l = n - 1;
        while (j > 0) {
            if (j >= i) {
                A[k][l] = true;
                j -= i;
            }
            i /= 2;
            l--;
        }
    }
}

void sum_vecb(unsigned int &s, vector<bool> &x) {
    // sums up the values of boolean x
    s = 0;
    for (auto &&i: x) {
        s += i;
    }
}

bool nonz_vec(vector<double> &x, double &prec) {
    return std::any_of(x.begin(), x.end(), [&prec](int i) { return i > prec || i < -prec; });
}

void vec_subtract(vector<double> &z, vector<double> &x, vector<double> &y) {
    // subtracts vector (double) y from vector (double) x
    for (unsigned int i = 0; i != x.size(); i++)
        z[i] = x[i] - y[i];
}

void sc_vec_prod(vector<double> &y, double a, vector<double> &x) {
    for (unsigned int i = 0; i < x.size(); i++)
        y[i] = a * x[i];
}