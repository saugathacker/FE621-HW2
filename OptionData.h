#pragma once

#include <string>
#include <chrono>

// struct to hold all the information of the option chain and the calculated IVs and Greeks
struct OptionData
{
  // information from the downloaded data
  std::string expiration;
  double timeToMaturity;
  double strike;
  std::string optionType;
  double lastPrice;
  double bid;
  double ask;
  double impliedVolatility;
  bool inTheMoney;

  // calculated Implied Vol
  double bisectionImpliedVol;
  double bisectionTime;

  // calculated parity and bs price
  double bs_price;
  double binom_price;
  double american_binom_price;
  double trinom_price;
  double american_trinom_price;

  // Constructor for initialization
  OptionData(const std::string &exp, double ttm, double strk, const std::string &type,
             double lp, double b, double a, double iv, bool itm)
      : expiration(exp), timeToMaturity(ttm), strike(strk), optionType(type),
        lastPrice(lp), bid(b), ask(a),
        impliedVolatility(iv), inTheMoney(itm), bisectionImpliedVol(0), bisectionTime(0) {}

  void calculate_iv_and_greeks(double spotPrice, double interestRate);
  void calculate_bs_price(double spot, double rate);
  void calculate_binom_tree_price(double spot, double rate, int steps);
  void calculate_american_binom_tree_price(double spot, double rate, int steps);
};