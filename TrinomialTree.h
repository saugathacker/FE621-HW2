#pragma once
#include "BinomialTree.h"

class TrinomialTree : public BinomialTree
{
public:
    TrinomialTree(double spot, double strike, double time_to_maturity,
                  double interest_rate, double sigma, OptionType option_type)
        : BinomialTree(spot, strike, time_to_maturity, interest_rate, sigma, option_type) {};
    double operator()(int steps) const override;
};