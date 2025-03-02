#pragma once
#include "BinomialTree.h"

// inherits the Binomial tree and only override the () operator to price option
class TrinomialTree : public BinomialTree
{
public:
    TrinomialTree(double spot, double strike, double time_to_maturity,
                  double interest_rate, double sigma, OptionType option_type)
        : BinomialTree(spot, strike, time_to_maturity, interest_rate, sigma, option_type) {};
    double operator()(int steps) const override;
    double get_barrier_option_price(int steps, BarrierOptionType barrier_type, double barrier_level);
};