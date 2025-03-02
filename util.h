#pragma once
#include <map>
#include <set>
#include "BSM.h"
#include "BinomialTree.h"

// Normal CDF function (declaration)
double norm_cdf(double x);

// Normal PDF function
double norm_pdf(double x);

std::map<int, double> absolute_error(BSM &bs, BinomialTree &bt);

// Bisection Method
double bisection_method(BSM &bs, double market_price, bool debug = false);

std::set<double> early_exercise_strikes(double q, bool is_cont);