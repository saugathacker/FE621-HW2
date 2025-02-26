#pragma once
#include <iostream>
#include "OptionType.h"

class BinomialTree
{
public:
    BinomialTree(double spot, double strike, double time_to_maturity,
                 double interest_rate, double sigma, OptionType option_type);
    virtual double operator()(int steps) const;
    virtual double payoff(double spot) const;
    double get_spot() const { return spot_; }
    double get_strike() const { return strike_; }
    double get_time_to_maturity() const { return time_to_maturity_; }
    double get_interest_rate() const { return interest_rate_; }
    double get_sigma_() const { return sigma_; }
    OptionType get_option_type() const { return option_type_; }
    friend std::ostream &operator<<(std::ostream &os, const BinomialTree &bt);

protected:
    double spot_, strike_;
    double time_to_maturity_, interest_rate_, sigma_;
    OptionType option_type_;
};
