#include "TrinomialTree.h"

#include <iostream> // Required for debugging

double TrinomialTree::operator()(int steps) const
{
    double dt = time_to_maturity_ / steps;
    double drift, dx, pu, pm, pd;

    drift = (interest_rate_ - 0.5 * sigma_ * sigma_);
    dx = sigma_ * std::sqrt(3 * dt);

    pu = 0.5 * ((sigma_ * sigma_ * dt + drift * drift * dt * dt) / (dx * dx) + (drift * dt) / dx);
    pd = 0.5 * ((sigma_ * sigma_ * dt + drift * drift * dt * dt) / (dx * dx) - (drift * dt) / dx);
    pm = 1 - pu - pd;
    double discount_factor = std::exp(-interest_rate_ * dt);

    std::vector<double> optionPrices(2 * steps + 1, 0.0);

    // Compute terminal node payoffs
    for (int i = 0; i <= 2 * steps; i++)
    {
        double log_spot_node = std::log(spot_) + (steps - i) * dx;
        double spot_node = std::exp(log_spot_node);
        optionPrices[i] = payoff(spot_node);
    }

    // Backward induction with trinomial probabilities
    for (int j = steps - 1; j >= 0; j--)
    {
        for (int i = 0; i <= 2 * j; i++)
        {
            optionPrices[i] = discount_factor * (pu * optionPrices[i] + pm * optionPrices[i + 1] + pd * optionPrices[i + 2]);

            // Early exercise check for American options
            if (option_type_ == OptionType::AmericanCall || option_type_ == OptionType::AmericanPut)
            {
                double early_exercise = payoff(std::exp(std::log(spot_) + (j - i) * dx));
                optionPrices[i] = std::max(optionPrices[i], early_exercise);
            }
        }
    }
    return optionPrices[0];
}
