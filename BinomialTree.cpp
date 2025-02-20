#include "BinomialTree.h"
#include <cmath>
#include <vector>

BinomialTree::BinomialTree(double spot, double strike, double time_to_maturity,
                           double interest_rate, double sigma, OptionType option_type)
    : spot_(spot), strike_(strike), time_to_maturity_(time_to_maturity),
      interest_rate_(interest_rate), sigma_(sigma), option_type_(option_type) {}

double BinomialTree::operator()(int steps) const
{
    double dt = time_to_maturity_ / steps;
    double dx, pu, pd;

    dx = std::sqrt(std::pow((interest_rate_ - 0.5 * sigma_ * sigma_), 2) * dt * dt + sigma_ * sigma_ * dt);
    pu = 0.5 + ((interest_rate_ - 0.5 * sigma_ * sigma_) * dt) / (2.0 * dx);
    pd = 1 - pu;
    double discount_factor = std::exp(-interest_rate_ * dt); // dicount factor for each step
    std::vector<double> optionPrices(steps + 1);

    for (int i = 0; i <= steps; i++)
    {
        double log_spot_node = std::log(spot_) + (steps - 2 * i) * dx;
        double spot_node = std::exp(log_spot_node);
        optionPrices[i] = payoff(spot_node);
    }

    for (int j = steps - 1; j >= 0; j--)
    {
        for (int i = 0; i <= j; i++)
        {
            optionPrices[i] = discount_factor * (pu * optionPrices[i] + pd * optionPrices[i + 1]);
        }
    }

    return optionPrices[0];
}

double BinomialTree::payoff(double spot) const
{
    int phi = static_cast<int>(option_type_);

    return std::max(spot * phi - phi * strike_, 0.0);
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
       << "Option Type: " << (bt.option_type_ == OptionType::Call ? "Call" : "Put") << "\n";
    return os;
}
