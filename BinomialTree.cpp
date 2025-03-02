#include "BinomialTree.h"
#include <cmath>
#include <vector>

BinomialTree::BinomialTree(double spot, double strike, double time_to_maturity,
                           double interest_rate, double sigma, OptionType option_type)
    : spot_(spot), strike_(strike), time_to_maturity_(time_to_maturity),
      interest_rate_(interest_rate), sigma_(sigma), option_type_(option_type) {}

// overloading the () operator to calculate the option price by taking steps as input
double BinomialTree::operator()(int steps) const
{
    double dt = time_to_maturity_ / steps;
    double drift, dx, pu, pd;

    drift = (interest_rate_ - 0.5 * sigma_ * sigma_);
    // implent the formula dor dx, pu and pd
    dx = std::sqrt(std::pow(drift, 2) * dt * dt + sigma_ * sigma_ * dt);
    pu = 0.5 + (drift * dt) / (2.0 * dx);
    pd = 1 - pu;
    double discount_factor = std::exp(-interest_rate_ * dt); // dicount factor for each step
    std::vector<double> optionPrices(steps + 1);

    for (int i = 0; i <= steps; i++)
    {
        // payoff at the terminal node
        double log_spot_node = std::log(spot_) + (steps - 2 * i) * dx;
        double spot_node = std::exp(log_spot_node);
        optionPrices[i] = payoff(spot_node);
    }

    // backward induction
    for (int j = steps - 1; j >= 0; j--)
    {
        for (int i = 0; i <= j; i++)
        {

            optionPrices[i] = discount_factor * (pu * optionPrices[i] + pd * optionPrices[i + 1]);

            // check early exercise for American options
            if (option_type_ == OptionType::AmericanCall || option_type_ == OptionType::AmericanPut)
            {
                double early_exercise = payoff(std::exp(std::log(spot_) + (j - 2 * i) * dx));
                optionPrices[i] = std::max(optionPrices[i], early_exercise);
            }
        }
    }

    return optionPrices[0];
}

double BinomialTree::payoff(double spot) const
{
    if (option_type_ == OptionType::EuropeanCall || option_type_ == OptionType::AmericanCall)
        return std::max(spot - strike_, 0.0);
    else
        return std::max(strike_ - spot, 0.0);
}

// Overloaded output stream operator
std::ostream &operator<<(std::ostream &os, const BinomialTree &bt)
{
    os << "Binomial Tree Parameters:\n"
       << "Spot Price: " << bt.spot_ << "\n"
       << "Strike Price: " << bt.strike_ << "\n"
       << "Time to Maturity: " << bt.time_to_maturity_ << " years\n"
       << "Interest Rate: " << bt.interest_rate_ * 100 << "%\n"
       << "Volatility: " << bt.sigma_ * 100 << "%\n"
       << "Option Type: " << (bt.option_type_ == OptionType::AmericanCall || bt.option_type_ == OptionType::EuropeanCall ? "Call" : "Put") << "\n";
    return os;
}
