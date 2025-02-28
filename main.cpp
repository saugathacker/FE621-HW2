#include "BinomialTree.h"
#include "TrinomialTree.h"
#include "BarrierBSM.h"

int main()
{

    // Define option parameters
    double spot = 100.0;         // Current stock price
    double strike = 100.0;       // Strike price
    double expiry = 1.0;         // Time to maturity (in years)
    double interest_rate = 0.05; // Risk-free interest rate
    double volatility = 0.4;     // Volatility (sigma)
    int steps = 10000;           // Number of steps in binomial tree

    // // Create BinomialTree objects for Call and Put options
    // BinomialTree euro_call_tree(spot, strike, expiry, interest_rate, volatility, OptionType::EuropeanCall);
    // BinomialTree euro_put_tree(spot, strike, expiry, interest_rate, volatility, OptionType::EuropeanPut);
    // BinomialTree ameri_call_tree(spot, strike, expiry, interest_rate, volatility, OptionType::AmericanCall);
    // BinomialTree ameri_put_tree(spot, strike, expiry, interest_rate, volatility, OptionType::AmericanPut);

    // // Compute option prices using BinomialTree
    // double call_price = euro_call_tree(steps);
    // double put_price = euro_put_tree(steps);

    // // Output results
    // std::cout << "European Call Option Price: " << call_price << std::endl;
    // std::cout << "European Put Option Price: " << put_price << std::endl;

    // call_price = ameri_call_tree(steps);
    // put_price = ameri_put_tree(steps);

    // // Output results
    // std::cout << "American Call Option Price: " << call_price << std::endl;
    // std::cout << "American Put Option Price: " << put_price << std::endl;

    // std::cout << "===================================================" << std::endl;
    // std::cout << "Trinomial Tree" << std::endl;

    // // Create Trinomial Tree objects for Call and Put options
    // TrinomialTree trinom_euro_call_tree(spot, strike, expiry, interest_rate, volatility, OptionType::EuropeanCall);
    // TrinomialTree trinom_euro_put_tree(spot, strike, expiry, interest_rate, volatility, OptionType::EuropeanPut);
    // TrinomialTree trinom_ameri_call_tree(spot, strike, expiry, interest_rate, volatility, OptionType::AmericanCall);
    // TrinomialTree trinom_ameri_put_tree(spot, strike, expiry, interest_rate, volatility, OptionType::AmericanPut);

    // // Compute option prices using Trinomial Tree
    // double call_price3 = trinom_euro_call_tree(steps);
    // double put_price3 = trinom_euro_put_tree(steps);

    // // Output results
    // std::cout << "European Call Option Price: " << call_price3 << std::endl;
    // std::cout << "European Put Option Price: " << put_price3 << std::endl;

    // call_price3 = trinom_ameri_call_tree(steps);
    // put_price3 = trinom_ameri_put_tree(steps);

    // // Output results
    // std::cout << "American Call Option Price: " << call_price3 << std::endl;
    // std::cout << "American Put Option Price: " << put_price3 << std::endl;

    // Construct a trinomial tree to calculate the price of an European Up-and-Out
    // call option. Use S0 = 10, strike K = 10, maturity T = 0.3, volatility σ = 0.2,
    // short rate r = 0.01, dividends δ = 0 and barrier H = 11

    double s0 = 10;
    double k = 10;
    double maturity = 0.3;
    double sigma = 0.2;
    double r = 0.01;

    double barrier_level = 11;

    TrinomialTree trinom_euro_call_tree2(s0, k, maturity, r, sigma, OptionType::EuropeanCall);

    double trinom_up_and_in_price = trinom_euro_call_tree2.get_barrier_option_price(steps, BarrierOptionType::UpAndIn, barrier_level);
    double trinom_up_and_out_price = trinom_euro_call_tree2.get_barrier_option_price(steps, BarrierOptionType::UpAndOut, barrier_level);

    std::cout << "Trinomial Tree Prices" << std::endl;
    std::cout << "Up-and-In Barrier Call Price: " << trinom_up_and_in_price << std::endl;
    std::cout << "Up-and-Out Barrier Call Price: " << trinom_up_and_out_price << std::endl;

    BarrierBSM barrier_model(s0, k, r, maturity);

    double up_and_in_price = barrier_model.barrier_call_price(sigma, BarrierOptionType::UpAndIn, barrier_level);
    double up_and_out_price = barrier_model.barrier_call_price(sigma, BarrierOptionType::UpAndOut, barrier_level);

    std::cout << "===================================================" << std::endl;
    std::cout << "Barrier BSM Prices" << std::endl;
    std::cout << "Up-and-In Barrier Call Price: " << up_and_in_price << std::endl;
    std::cout << "Up-and-Out Barrier Call Price: " << up_and_out_price << std::endl;

    return 0;
}