#include "util.h"
#include <cmath>
#include <numbers>
#include <array>

// Implementation of the Normal CDF function
double norm_cdf(double x)
{
    return (1.0 + std::erf(x / std::numbers::sqrt2)) / 2;
}

// Implementation of the Normal PDF function
double norm_pdf(double x)
{
    return std::exp(-0.5 * x * x) / std::sqrt(2 * std::numbers::pi);
}

std::map<int, double> absolute_error(BSM &bs, BinomialTree &bt)
{
    std::array<int, 12> steps_list = {10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 350, 400};
    std::map<int, double> result;
    double vol = 0.4;
    // iterate all the choices for steps
    for (int steps : steps_list)
    {
        double error = std::fabs(bs(vol) - bt(steps));
        result[steps] = error;
    }

    return result;
}

double bisection_method(BSM &bs, double market_price, bool debug)
{
    double a = 0.0001;
    double b = 3.0;
    double c = (a + b) / 2;
    double epsilon = 1e-06;
    double tol = bs(c) - market_price;
    int max_iter = 1000;
    int iter = 0;

    // **Check if the root is even bracketed**
    double fa = bs(a) - market_price;
    double fb = bs(b) - market_price;

    if (debug)
        std::cout << "bs(a) " << bs(a) << "bs(b) " << bs(b) << std::endl;

    if (fa * fb > 0)
    {
        if (debug)
            std::cout << "Warning: Root is not in range! Returning best estimate.\n";
        return std::min(std::max(epsilon, a), b);
    }

    if (debug)
        std::cout << "Starting Bisection Method...\n";

    while (std::abs(tol) > epsilon && iter < max_iter)
    {
        tol = bs(c) - market_price;

        if (tol < 0)
        {
            a = c;
        }
        else
        {
            b = c;
        }

        c = (a + b) / 2;
        iter++;

        // **Optional debug logging**
        if (debug && iter % 10 == 0)
        {
            std::cout << "Iteration " << iter << ": IV estimate = " << c << " | Tolerance = " << tol << std::endl;
        }
    }

    if (debug)
        std::cout << "Ending Bisection Method after " << iter << " iterations.\n";

    // **Ensure valid output within IV range**
    c = std::max(epsilon, std::min(c, b));

    return c;
}

std::set<double> early_exercise_strikes(double q, bool is_cont)
{
    double S0 = 40;
    double N = 2;
    double T = 0.5;
    double u = 1.2;
    double d = 0.9;
    double r = 0.04;
    double dt = T / N;
    double discount = std::exp(-r * dt);
    double p = is_cont ? (std::exp((r - q) * dt) - d) / (u - d) : (std::exp(r * dt) - d) / (u - d);
    double D = S0 * q;

    std::set<double> early_exercise_strikes;

    // generate possible strike prices in a reasonable range
    for (double K = 1; K <= 100.0; K++)
    {
        // build the binomial tree
        std::vector<std::vector<double>> stock_tree(N + 1, std::vector<double>(N + 1, 0.0));
        std::vector<std::vector<double>> option_tree(N + 1, std::vector<double>(N + 1, 0.0));

        if (is_cont)
        {
            for (int j = 0; j <= N; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    stock_tree[i][j] = S0 * std::pow(u, j - i) * std::pow(d, i);
                }
            }
        }
        else
        {
            // time 0
            stock_tree[0][0] = S0;

            // time 1
            stock_tree[1][0] = (S0 * u) - D;
            stock_tree[1][1] = (S0 * d) - D;

            // time 2
            stock_tree[2][0] = stock_tree[1][0] * u;
            stock_tree[2][1] = stock_tree[1][0] * d;
            stock_tree[2][2] = stock_tree[1][1] * d;
        }

        // compute option prices at expiration
        for (int i = 0; i <= N; ++i)
        {
            option_tree[i][N] = std::max(stock_tree[i][N] - K, 0.0);
        }

        // backward induction to check early exercise
        bool early_exercise_found = false;
        for (int j = N - 1; j >= 0; --j)
        {
            for (int i = 0; i <= j; ++i)
            {
                double intrinsic_value = std::max(stock_tree[i][j] - K, 0.0);
                double expected_value = discount * (p * option_tree[i][j + 1] + (1 - p) * option_tree[i + 1][j + 1]);

                option_tree[i][j] = std::max(intrinsic_value, expected_value);

                if (intrinsic_value > expected_value)
                {
                    early_exercise_found = true;
                }
            }
        }

        if (early_exercise_found)
        {
            early_exercise_strikes.insert(K);
        }
    }

    return early_exercise_strikes;
}