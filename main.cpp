#include "BinomialTree.cpp"

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
    BinomialTree call_tree(spot, strike, expiry, interest_rate, volatility, OptionType::Call);
    BinomialTree put_tree(spot, strike, expiry, interest_rate, volatility, OptionType::Put);

    // Compute option prices using BinomialTree
    double call_price = call_tree(steps);
    double put_price = put_tree(steps);

    // Output results
    std::cout << "European Call Option Price: " << call_price << std::endl;
    std::cout << "European Put Option Price: " << put_price << std::endl;

    return 0;
}