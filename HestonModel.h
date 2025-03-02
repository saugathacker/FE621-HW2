#pragma once

#include <vector>

#include <cmath>
#include <memory>
#include <iostream>
#include "OptionType.h"
#include <iomanip>

using namespace std;

// Variance node to store variance, transformed variance and probs
struct VarianceNode
{
    double xt;
    double vt;
    vector<double> probs;

    VarianceNode(double x, double v) : xt(x), vt(v), probs(3, 0.0) {}
};

// stock node to store stock, transformed stock and probs
struct StockNode
{
    double yt;
    double st;
    vector<double> probs;

    StockNode(double y, double s) : yt(y), st(s), probs(3, 0.0) {}
};

class HestonModel
{
public:
    HestonModel(int steps, double S0, double V0, double kappa, double theta,
                double sigma, double rho, double r, double T, double K, OptionType option_type);
    double price_option();
    void display_vol_tree();
    void display_stock_tree();
    void construct_vol_tree();
    void construct_stock_tree();

private:
    int steps;
    double S0, V0, kappa, theta, sigma, rho, r, T, K;
    double dt, b, threshold;
    OptionType option_type;
    vector<vector<shared_ptr<VarianceNode>>> vol_tree;
    vector<vector<shared_ptr<StockNode>>> stock_tree;
    void compute_b();
    vector<double> compute_vol_prob(double xt);
    vector<double> compute_stock_prob(double Vt, double Yt);
    double payoff(double s);
};
