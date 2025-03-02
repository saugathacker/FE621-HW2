#include "Ticker.h"
#include <fstream>
#include <iostream>

// Constructor
Ticker::Ticker(const std::string &name, double spot, double rate)
    : tickerName(name), spotPrice(spot), interestRate(rate) {}

// add new option data to the existing Ticker object
void Ticker::addOptionData(std::unique_ptr<OptionData> option)
{
    options.push_back(std::move(option));
}

// find an option based on strike, expiration, and type
OptionData *Ticker::findOption(double strike, const std::string &expiration, const std::string &optionType) const
{
    for (const auto &option : options)
    {
        if (option->strike == strike && option->expiration == expiration && option->optionType == optionType)
        {
            return option.get(); // return raw pointer (safe since unique_ptr manages memory)
        }
    }
    return nullptr; // return nullptr if no matching option is found
}

void Ticker::calculate_implied_vols_and_bs_price()
{
    for (auto &option : options)
    {
        option->calculate_iv_and_greeks(spotPrice, interestRate); // each option calculates its IV
        option->calculate_bs_price(spotPrice, interestRate);
    }
}

void Ticker::calculate_binom_tree_price()
{
    for (auto &option : options)
    {
        option->calculate_binom_tree_price(spotPrice, interestRate, 1000);
        option->calculate_american_binom_tree_price(spotPrice, interestRate, 1000);
    }
}

// calcullate the iv using tree methods
void Ticker::calculate_tree_iv()
{
    for (auto &option : options)
    {
        option->calculate_binom_iv(spotPrice, interestRate);
        option->calculate_trinom_iv(spotPrice, interestRate);
    }
}

// implentaion of write to csv all the option data (observed and calculated) of this Ticker
void Ticker::write_to_csv(const std::string &filename) const
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    // **Write CSV Header**
    file << "Ticker,Expiration,TimeToMaturity,Strike,OptionType,LastPrice,"
         << "Bid,Ask,ImpliedVolatility,BisectionIV,BisectionTime,Bs_price,Binom_price,American_binom_price,Binom_bisection_iv, Trinom_bisection_iv,InTheMoney\n";

    // **Write Option Data**
    for (const auto &option : options)
    {
        file << tickerName << "," // Add ticker symbol
             << option->expiration << ","
             << option->timeToMaturity << ","
             << option->strike << ","
             << option->optionType << ","
             << option->lastPrice << ","
             << option->bid << ","
             << option->ask << ","
             << option->impliedVolatility << ","
             << option->bisectionImpliedVol << ","
             << option->bisectionTime << ","
             << option->bs_price << ","
             << option->binom_price << ","
             << option->american_binom_price << ","
             << option->binomBisectioIV << ","
             << option->trinomBisectionIV << ","
             << (option->inTheMoney ? "True" : "False") << "\n";
    }

    file.close();
    std::cout << "CSV file written successfully: " << filename << std::endl;
}