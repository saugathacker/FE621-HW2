#include "BinomialTree.h"
#include "TrinomialTree.h"

int main()
{

    // Define option parameters
    double spot = 100.0;         // Current stock price
    double strike = 100.0;       // Strike price
    double expiry = 1.0;         // Time to maturity (in years)
    double interest_rate = 0.05; // Risk-free interest rate
    double volatility = 0.2;     // Volatility (sigma)
    int steps = 100;             // Number of steps in binomial tree

    // Create BinomialTree objects for Call and Put options
    BinomialTree euro_call_tree(spot, strike, expiry, interest_rate, volatility, OptionType::EuropeanCall);
    BinomialTree euro_put_tree(spot, strike, expiry, interest_rate, volatility, OptionType::EuropeanPut);
    BinomialTree ameri_call_tree(spot, strike, expiry, interest_rate, volatility, OptionType::AmericanCall);
    BinomialTree ameri_put_tree(spot, strike, expiry, interest_rate, volatility, OptionType::AmericanPut);

    // Compute option prices using BinomialTree
    double call_price = euro_call_tree(steps);
    double put_price = euro_put_tree(steps);

    // Output results
    std::cout << "European Call Option Price: " << call_price << std::endl;
    std::cout << "European Put Option Price: " << put_price << std::endl;

    call_price = ameri_call_tree(steps);
    put_price = ameri_put_tree(steps);

    // Output results
    std::cout << "American Call Option Price: " << call_price << std::endl;
    std::cout << "American Put Option Price: " << put_price << std::endl;

    std::cout << "===================================================" << std::endl;
    std::cout << "Trinomial Tree" << std::endl;

    // Create Trinomial Tree objects for Call and Put options
    TrinomialTree trinom_euro_call_tree(spot, strike, expiry, interest_rate, volatility, OptionType::EuropeanCall);
    TrinomialTree trinom_euro_put_tree(spot, strike, expiry, interest_rate, volatility, OptionType::EuropeanPut);
    TrinomialTree trinom_ameri_call_tree(spot, strike, expiry, interest_rate, volatility, OptionType::AmericanCall);
    TrinomialTree trinom_ameri_put_tree(spot, strike, expiry, interest_rate, volatility, OptionType::AmericanPut);

    // Compute option prices using Trinomial Tree
    double call_price3 = trinom_euro_call_tree(steps);
    double put_price3 = trinom_euro_put_tree(steps);

    // Output results
    std::cout << "European Call Option Price: " << call_price3 << std::endl;
    std::cout << "European Put Option Price: " << put_price3 << std::endl;

    call_price3 = trinom_ameri_call_tree(steps);
    put_price3 = trinom_ameri_put_tree(steps);

    // Output results
    std::cout << "American Call Option Price: " << call_price3 << std::endl;
    std::cout << "American Put Option Price: " << put_price3 << std::endl;

    return 0;
}