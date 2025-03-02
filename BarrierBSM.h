#pragma once
#include "BSM.h"

// Class to implement explicit formula for barrier option using BlackScholes Model
class BarrierBSM
{
public:
    BarrierBSM(double spot, double strike, double interest_rate, double time_to_maturity)
        : spot_(spot), strike_(strike), interest_rate_(interest_rate), time_to_maturity_(time_to_maturity) {};
    double barrier_call_price(double vol, BarrierOptionType barrier_type, double barrier_level) const;

private:
    double spot_, strike_, interest_rate_, time_to_maturity_;
    double compute_d_BS(double A, double B, double vol) const;
    double black_scholes_price(double vol, OptionType option_type, double K, double S) const;
};
