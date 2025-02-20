#include "util.h"
#include <cmath>
#include <numbers>
#include <array>

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

std::unordered_map<int, double> absolute_error(BSM &bs, BinomialTree &bt)
{
    std::array<int, 12> steps_list = {10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 350, 400};
    std::unordered_map<int, double> result;
    double vol = 0.4;
    // iterate all the choices for steps
    for (int steps : steps_list)
    {
        double error = std::fabs(bs(vol) - bt(steps));
        result[steps] = error;
    }

    return result;
}