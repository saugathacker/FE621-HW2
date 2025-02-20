#pragma once
#include <unordered_map>
#include "BSM.h"
#include "BinomialTree.h"

// Normal CDF function (declaration)
double norm_cdf(double x);

// Normal PDF function
double norm_pdf(double x);

std::unordered_map<int, double> absolute_error(BSM &bs, BinomialTree &bt);