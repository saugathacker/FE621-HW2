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

double TrinomialTree::get_barrier_option_price(int steps, BarrierOptionType barrier_type, double barrier_level)
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
    double log_spot = std::log(spot_);
    // Compute terminal node payoffs
    for (int i = 0; i <= 2 * steps; i++)
    {
        double log_spot_node = log_spot + (steps - i) * dx;
        double spot_node = std::exp(log_spot_node);
        optionPrices[i] = payoff(spot_node);
    }

    // Backward induction with trinomial probabilities
    for (int j = steps - 1; j >= 0; j--)
    {
        for (int i = 0; i <= 2 * j; i++)
        {
            double current_stock_price = std::exp(log_spot + (j - i) * dx);
            double option_price = discount_factor * (pu * optionPrices[i] + pm * optionPrices[i + 1] + pd * optionPrices[i + 2]);

            bool is_activated = false;
            switch (barrier_type)
            {
            case BarrierOptionType::UpAndIn:
                is_activated = (current_stock_price >= barrier_level);
                break;
            case BarrierOptionType::UpAndOut:
                is_activated = (current_stock_price < barrier_level);
                break;
            case BarrierOptionType::DownAndIn:
                is_activated = (current_stock_price <= barrier_level);
                break;
            case BarrierOptionType::DownAndOut:
                is_activated = (current_stock_price > barrier_level);
                break;
            default:
                break;
            }

            if (!is_activated)
            {
                optionPrices[i] = 0.0; // Barrier breached (option knocked out)
            }
            else
            {
                if (option_type_ == OptionType::AmericanCall || option_type_ == OptionType::AmericanPut)
                {
                    double early_exercise = payoff(current_stock_price);
                    option_price = std::max(optionPrices[i], early_exercise);
                }
                optionPrices[i] = option_price;
            }
        }
    }
    return optionPrices[0];
}
