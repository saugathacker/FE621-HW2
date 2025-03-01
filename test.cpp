#include "HestonDraft.cpp"
#include "OptionType.h"

int main()
{
    int steps = 1000;
    double dt = 1 / 1000;
    double S0 = 100;
    double kappa = 1.5;
    double theta = 0.04;
    double V0 = 0.04;
    double sigma = 0.2;
    double rho = 0.8;
    double K = 100;
    double r = 0.01;

    HestonModel h1(steps, S0, V0, kappa, theta, sigma, rho, r, 1, K, OptionType::EuropeanCall);
}