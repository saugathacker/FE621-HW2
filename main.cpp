#include "BinomialTree.h"
#include "TrinomialTree.h"
#include "BarrierBSM.h"
#include <iomanip>

#include "FileHandler.h"
#include "Ticker.h"
#include "util.h"

int main()
{

    // import all the option data

    FileHandler fh(',', true);
    std::unordered_map<std::string, std::unique_ptr<Ticker>> tickers = fh.read_csv_into_ticker_object("NVDA_options_data.csv");
    for (const auto &[ticker, tickerObj] : tickers)
    {
        // calculate IV and binom tree price
        tickerObj->calculate_implied_vols_and_bs_price();
        tickerObj->calculate_binom_tree_price();
        tickerObj->write_to_csv("output_NVDA_test.csv");
    }

    // Define option parameters for Absolute Error
    double spot = 100.0;         // Current stock price
    double strike = 100.0;       // Strike price
    double expiry = 0.3;         // Time to maturity (in years)
    double interest_rate = 0.05; // Risk-free interest rate
    double volatility = 0.4;     // Volatility (sigma)
    int steps = 10000;           // Number of steps in binomial tree

    // Absolute Error
    // construct one BSM and one BinomTree Object
    BSM bs1(strike, spot, expiry, interest_rate, OptionType::EuropeanCall, 0.0);
    BinomialTree binom_tree1(spot, strike, expiry, interest_rate, volatility, OptionType::EuropeanCall);

    std::map<int, double> absolute_error_map = absolute_error(bs1, binom_tree1);

    std::cout << "================================\n";
    std::cout << "|   Steps   |   Absolute Error  |\n";
    std::cout << "================================\n";

    for (const auto &[steps, error] : absolute_error_map)
    {
        std::cout << "| " << std::setw(8) << steps << " | " << std::setw(16) << error << " |\n";
    }

    std::cout << "================================\n";

    // Construct a trinomial tree to calculate the price of an European Up-and-Out
    // call option. Use S0 = 10, strike K = 10, maturity T = 0.3, volatility σ = 0.2,
    // short rate r = 0.01, dividends δ = 0 and barrier H = 11

    double s0 = 10;
    double k = 10;
    double T = 0.3;
    double sigma = 0.2;
    double r = 0.01;
    double barrier_level = 11;

    TrinomialTree trinom_euro_call_tree2(s0, k, T, r, sigma, OptionType::EuropeanCall);

    double trinom_up_and_in_price = trinom_euro_call_tree2.get_barrier_option_price(10,
                                                                                    BarrierOptionType::UpAndIn, barrier_level);
    double trinom_up_and_out_price = trinom_euro_call_tree2.get_barrier_option_price(steps,
                                                                                     BarrierOptionType::UpAndOut, barrier_level);

    BarrierBSM barrier_model(s0, k, r, T);

    double up_and_in_price = barrier_model.barrier_call_price(sigma, BarrierOptionType::UpAndIn, barrier_level);
    double up_and_out_price = barrier_model.barrier_call_price(sigma, BarrierOptionType::UpAndOut, barrier_level);

    std::cout << "===============================================" << std::endl;
    std::cout << "|         Barrier Option Prices              |" << std::endl;
    std::cout << "===============================================" << std::endl;
    std::cout << "| " << std::setw(30) << "Method" << " | "
              << std::setw(20) << "Price" << " |" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "| " << std::setw(30) << "Trinomial Up-and-In Call" << " | "
              << std::setw(20) << trinom_up_and_in_price << " |" << std::endl;

    std::cout << "| " << std::setw(30) << "Trinomial Up-and-Out Call" << " | "
              << std::setw(20) << trinom_up_and_out_price << " |" << std::endl;

    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "| " << std::setw(30) << "BSM Up-and-In Call" << " | "
              << std::setw(20) << up_and_in_price << " |" << std::endl;

    std::cout << "| " << std::setw(30) << "BSM Up-and-Out Call" << " | "
              << std::setw(20) << up_and_out_price << " |" << std::endl;

    std::cout << "===============================================" << std::endl;

    // Trinomial Tree for American put option
    TrinomialTree american_put_trinom_tree(s0, k, T, T, sigma, OptionType::AmericanPut);

    double american_up_and_in_put = american_put_trinom_tree.get_barrier_option_price(steps, BarrierOptionType::UpAndIn, barrier_level);

    std::cout << "| " << std::setw(30) << "American Trinomial Up-and-In Put" << " | "
              << std::setw(20) << american_up_and_in_put << " |" << std::endl;

    // dividend assumptions

    double q_continuous = 0.02;
    double q_discrete = 0.03;

    std::set<double> strikes_continuous = early_exercise_strikes(q_continuous, true);
    std::set<double> strikes_discrete = early_exercise_strikes(q_discrete, false);

    std::cout << "=================================================" << std::endl;
    std::cout << "Early Exercise Strikes (Continuous 2% Dividend): ";
    for (double k : strikes_continuous)
    {
        std::cout << k << " ";
    }
    std::cout << std::endl;

    std::cout << "Early Exercise Strikes (Discrete 3% Dividend): ";
    for (double k : strikes_discrete)
    {
        std::cout << k << " ";
    }
    std::cout << std::endl;
    std::cout << "=================================================" << std::endl;

    // chooser option

    double S0_chooser = 100;
    double K_chooser = 100;
    double maturity_chooser = 0.3;
    double volatility_chooser = 0.15;
    double rate_chooser = 0.01;
    int steps_chooser = 1000;

    BinomialTree binom_chooser_euro(S0_chooser, K_chooser, maturity_chooser, rate_chooser, volatility_chooser, OptionType::EuropeanCall);
    BinomialTree binom_chooser_american(S0_chooser, K_chooser, maturity_chooser, rate_chooser, volatility_chooser, OptionType::AmericanCall);

    double euro_chooser = binom_chooser_euro.chooser_option(steps_chooser);
    double american_chooser = binom_chooser_american.chooser_option(steps_chooser);

    // Display parameters and results
    std::cout << std::fixed << std::setprecision(4); // Format output to 4 decimal places
    std::cout << "=============================================\n";
    std::cout << "          Chooser Option Pricing            \n";
    std::cout << "=============================================\n";
    std::cout << "Stock Price (S0_chooser):      " << S0_chooser << std::endl;
    std::cout << "Strike Price (K):      " << K_chooser << std::endl;
    std::cout << "Time to Maturity (T):  " << maturity_chooser << " years\n";
    std::cout << "Volatility (sigma):    " << volatility_chooser << std::endl;
    std::cout << "Risk-Free Rate (r):    " << rate_chooser << std::endl;
    std::cout << "Binomial Steps:        " << steps_chooser << std::endl;
    std::cout << "---------------------------------------------\n";
    std::cout << "European Chooser Option Price: " << euro_chooser << std::endl;
    std::cout << "American Chooser Option Price: " << american_chooser << std::endl;
    std::cout << "=============================================\n";

    std::unordered_map<std::string, std::unique_ptr<Ticker>> tickers_data1 = fh.read_csv_into_ticker_object2("options_data1.csv");
    std::unordered_map<std::string, std::unique_ptr<Ticker>> tickers_data2 = fh.read_csv_into_ticker_object2("options_data2.csv");

    // calculate implied vol for data1 and writing it into csv file
    for (const auto &[ticker, tickerObj] : tickers_data1)
    {
        // IV using BS
        tickerObj->calculate_implied_vols_and_bs_price();
        // IV using binom tree
        // IV using trinom tree
        tickerObj->calculate_tree_iv();
        std::string outputFileName = ticker + "_outputData1.csv";
        tickerObj->write_to_csv(outputFileName);

        // for each ticker calculate the option price using calculated implied volatitlity from previous day
        tickers_data2[ticker]->calculate_implied_vols_and_bs_price();
        tickers_data2[ticker]->calculate_tree_iv();
        outputFileName = ticker + "_outputData2.csv";
        tickers_data2[ticker]->write_to_csv(outputFileName);
    }
    return 0;
}