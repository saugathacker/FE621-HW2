#pragma once
#include <string>
#include <vector>
#include <memory>
#include "OptionData.h"

// class to store and manage options for a specific ticker
class Ticker
{
private:
    std::string tickerName; // stock ticker symbol
    double spotPrice;
    double interestRate;
    std::vector<std::unique_ptr<OptionData>> options; // unique ptr vector of OptionData

public:
    // constructor
    Ticker(const std::string &name, double spot, double rate);

    // add option data to the existing Ticker object
    void addOptionData(std::unique_ptr<OptionData> option);

    // getter for ticker name and option chain size
    std::string getTickerName() const { return tickerName; }
    double getOptionsSize() const { return options.size(); };
    // getter for spot price and interest rate
    double getSpotPrice() const { return spotPrice; };
    double getInterestRate() const { return interestRate; };

    // find and return a pointer to the OptionData that matches strike, expiration and type
    OptionData *findOption(double strike, const std::string &expiration, const std::string &optionType) const;

    // functions to calculate the implied vol, greeks, parity price and bs price
    void calculate_implied_vols_and_bs_price();
    void calculate_binom_tree_price();

    // function to write all options to a CSV file
    void write_to_csv(const std::string &filename) const;
};