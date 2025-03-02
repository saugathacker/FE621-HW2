#include "OptionData.h"
#include <iostream>
#include <cmath>
#include "util.h"
#include "BSM.h"
#include "BinomialTree.h"
#include "TrinomialTree.h"

// implementation of calculate iv and greeks
void OptionData::calculate_iv_and_greeks(double spotPrice, double interestRate)
{
    // Skip calculation if lastPrice, bid, and ask are all zero
    if ((lastPrice <= 0.0 || std::isnan(lastPrice)) &&
        (bid <= 0.0 || std::isnan(bid)) &&
        (ask <= 0.0 || std::isnan(ask)))
    {
        return; // Do not calculate IV if no valid market data exists
    }

    // Use mid-price if bid/ask exist, otherwise use last price
    double market_price = (bid > 0 && ask > 0) ? (bid + ask) / 2 : lastPrice;

    // Create Black-Scholes model instance
    BSM bs(strike, spotPrice, timeToMaturity, interestRate,
           (optionType == "Call" ? OptionType::EuropeanCall : OptionType::EuropeanPut));

    // std::cout << bs << std::endl;
    // std::cout << "Finding root with market price: " << market_price << std::endl;
    // Measure time for Bisection Method
    auto start_bisect = std::chrono::high_resolution_clock::now();
    bisectionImpliedVol = bisection_method(bs, market_price, false);
    auto end_bisect = std::chrono::high_resolution_clock::now();
    bisectionTime = std::chrono::duration<double, std::milli>(end_bisect - start_bisect).count();
}

void OptionData::calculate_bs_price(double spot, double rate)
{
    BSM bs_model(strike, spot, timeToMaturity, rate,
                 (optionType == "Call" ? OptionType::EuropeanCall : OptionType::EuropeanPut), 0.0);

    bs_price = bs_model(bisectionImpliedVol);
}

void OptionData::calculate_binom_tree_price(double spot, double rate, int steps)
{
    BinomialTree binom_tree(spot, strike, timeToMaturity, rate, bisectionImpliedVol, (optionType == "Call" ? OptionType::EuropeanCall : OptionType::EuropeanPut));
    binom_price = binom_tree(steps);
}

void OptionData::calculate_american_binom_tree_price(double spot, double rate, int steps)
{
    BinomialTree binom_tree(spot, strike, timeToMaturity, rate, bisectionImpliedVol, (optionType == "Call" ? OptionType::AmericanCall : OptionType::AmericanPut));
    american_binom_price = binom_tree(steps);
}

// calculating IV using tree methods
void OptionData::calculate_binom_iv(double spotPrice, double interestRate)
{
    // Skip calculation if lastPrice, bid, and ask are all zero
    if ((lastPrice <= 0.0 || std::isnan(lastPrice)) &&
        (bid <= 0.0 || std::isnan(bid)) &&
        (ask <= 0.0 || std::isnan(ask)))
    {
        return; // Do not calculate IV if no valid market data exists
    }

    // Use mid-price if bid/ask exist, otherwise use last price
    double market_price = (bid > 0 && ask > 0) ? (bid + ask) / 2 : lastPrice;

    BinomialTree bt(spotPrice, strike, timeToMaturity, interestRate, 1.5,
                    (optionType == "Call" ? OptionType::EuropeanCall : OptionType::EuropeanPut));

    binomBisectioIV = tree_bisection(bt, market_price);
}

void OptionData::calculate_trinom_iv(double spotPrice, double interestRate)
{
    // Skip calculation if lastPrice, bid, and ask are all zero
    if ((lastPrice <= 0.0 || std::isnan(lastPrice)) &&
        (bid <= 0.0 || std::isnan(bid)) &&
        (ask <= 0.0 || std::isnan(ask)))
    {
        return; // Do not calculate IV if no valid market data exists
    }

    // Use mid-price if bid/ask exist, otherwise use last price
    double market_price = (bid > 0 && ask > 0) ? (bid + ask) / 2 : lastPrice;

    TrinomialTree trinom_tree(spotPrice, strike, timeToMaturity, interestRate, 1.5,
                              (optionType == "Call" ? OptionType::EuropeanCall : OptionType::EuropeanPut));

    trinomBisectionIV = tree_bisection(trinom_tree, market_price);
}