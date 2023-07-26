//
// Created by secci on 25/07/2023.
//
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "nucleolus.h"
#include "nucleolus_test.h"

TEST_CASE("Simple nucleolus of 3 players") {
    unsigned short int n = 3;
    double x[] = {0, 0, 0};
    double v[] = {0, 0, 0, 5, 5, 1, 9};
    primal_dual(v, x, n);
    double x_0 = round(x[0] * 100.0) / 100.0;
    CHECK(x_0 == 1.33);
}
