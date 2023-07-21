#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <time.h>
#include <fstream>
#include <glpk.h>

using namespace std;

void swap_ith_and_firstone(vector<vector<double>> &rref, vector<int> &ones, vector<int> &nonz, unsigned int &i);

void swap_ith_and_firstnz(vector<vector<double>> &rref, vector<int> &nonz, unsigned int &i);

void rowechform_piv2(vector<vector<double>> &rref, unsigned int &i, unsigned short int &n);

void
rowechform_piv(vector<vector<double>> &rref, vector<int> &nonz, unsigned int &i, unsigned short int &j, unsigned int &k,
               unsigned short int &n);

void rowechform_loop(vector<vector<double>> &rref, vector<bool> &J, unsigned int &i, unsigned short int &j,
                     unsigned short int &rank, double &prec, unsigned short int &n);

void rowechform(vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &B, unsigned short int &n,
                unsigned short int &rank);

bool binrank(vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &b, unsigned short int &n);

void A_mx(vector<vector<bool>> &A, unsigned short int &n, unsigned int &s);

void de2bi(unsigned int &k, vector<bool> &a, unsigned short int &n);

void sum_vecb(unsigned int &s, vector<bool> &x);

bool nonz_vec(vector<double> &x, double &prec);

void sc_vec_prod(vector<double> &y, double a, vector<double> &x);

void vec_subtract(vector<double> &z, vector<double> &x, vector<double> &y);

double cpuTime();