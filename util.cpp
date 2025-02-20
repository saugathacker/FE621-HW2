#include "util.h"
#include <cmath>
#include <numbers>

// Implementation of the Normal CDF function
double norm_cdf(double x)
{
    return (1.0 + std::erf(x / std::numbers::sqrt2)) / 2;
}

// Implementation of the Normal PDF function
double norm_pdf(double x)
{
    return std::exp(-0.5 * x * x) / std::sqrt(2 * std::numbers::pi);
}