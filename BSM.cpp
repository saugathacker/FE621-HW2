#include "BSM.h"
#include "util.h"

// Constructor Implementation
BSM::BSM(double strike, double spot, double time_to_maturity,
         double interest_rate, OptionType option_type, double dividend_yield)
    : strike_(strike), spot_(spot), time_to_maturity_(time_to_maturity),
      interest_rate_(interest_rate), dividend_yield_(dividend_yield), option_type_(option_type) {}

// Functor Implementation which takes volatility as an input and outputs the option price
double BSM::operator()(double vol) const
{

    using std::exp;
    // getting d1 and d2 through helper function
    auto norm_args = compute_norm_args_(vol);
    double d1 = norm_args[0];
    double d2 = norm_args[1];

    // calculating phi to help determine optiontype
    int phi = (option_type_ == OptionType::EuropeanCall || option_type_ == OptionType::AmericanCall) ? 1 : -1;

    double nD1 = norm_cdf(phi * d1);
    double nD2 = norm_cdf(phi * d2);

    double discountFactor = exp(-interest_rate_ * time_to_maturity_);
    double dividendFactor = exp(-dividend_yield_ * time_to_maturity_);

    // call option payoff when phi is 1 and put option payoff when phi is -1
    double payoff = phi * (spot_ * dividendFactor * nD1 - strike_ * discountFactor * nD2);

    return payoff;
}

// Helper function to compute d1 and d2
std::array<double, 2> BSM::compute_norm_args_(double vol) const
{

    // implementing the formula for d1 and d2
    double log_term = std::log(spot_ / strike_);

    double drift_term = (interest_rate_ - dividend_yield_ + 0.5 * vol * vol) * time_to_maturity_;

    double denom = vol * std::sqrt(time_to_maturity_);

    double d1 = (log_term + drift_term) / denom;

    double d2 = d1 - denom;

    return {d1, d2}; // Explicitly constructing std::array
}