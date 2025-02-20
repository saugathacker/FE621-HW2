#pragma once
#include <array>
#include <iostream>
#include "OptionType.h"

// Blackscholes class
class BSM
{
public:
    BSM(double strike, double spot, double time_to_maturity,
        double interest_rate, OptionType option_type, double dividend_yield = 0);
    double operator()(double vol) const;

    double get_spot() const { return spot_; }
    double get_strike() const { return strike_; }
    double get_time_to_maturity() const { return time_to_maturity_; }
    double get_interest_rate() const { return interest_rate_; }
    OptionType get_option_type() const { return option_type_; }
    friend std::ostream &operator<<(std::ostream &os, const BSM &bs);

private:
    // member variables spot price, strike price, time to maturity,
    // interest rate, dividend yield, option type
    double spot_, strike_, time_to_maturity_;
    double interest_rate_, dividend_yield_;
    OptionType option_type_;
    std::array<double, 2> compute_norm_args_(double vol) const;
};