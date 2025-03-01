#pragma once

#include <vector>
#include <cmath>
#include <memory>
#include <iostream>
#include "OptionType.h"

using namespace std;

struct BivariateNode
{
    double S;    // Stock Price
    double V;    // Variance
    double prob; // Joint Probability

    BivariateNode(double S, double V, double p) : S(S), V(V), prob(p) {}
};

struct VarianceNode
{
    double xt;
    double vt;
    vector<double> probs;

    VarianceNode(double x, double v) : xt(x), vt(v), probs(3, 0.0) {}
};

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
    void construct_bivariate_tree();
    vector<vector<shared_ptr<BivariateNode>>> bivariate_tree;

private:
    int steps;
    double S0, V0, kappa, theta, sigma, rho, r, T, K;
    double dt, b, threshold;
    OptionType option_type;
    vector<vector<shared_ptr<VarianceNode>>> vol_tree;
    vector<vector<shared_ptr<StockNode>>> stock_tree;
    void compute_b();
    void construct_vol_tree();
    vector<double> compute_vol_prob(double xt);
    vector<double> compute_stock_prob(double Vt, double Yt);
    void construct_stock_tree();
};

HestonModel::HestonModel(int steps, double S0, double V0, double kappa, double theta,
                         double sigma, double rho, double r, double T, double K,
                         OptionType option_type) : steps(steps), S0(S0), V0(V0),
                                                   kappa(kappa), theta(theta),
                                                   sigma(sigma), rho(rho), r(r),
                                                   T(T), K(K), option_type(option_type)
{
    dt = T / steps;
    threshold = 1e-06;
    compute_b();
    construct_vol_tree();
    construct_stock_tree();
    construct_bivariate_tree();
}

void HestonModel::compute_b()
{
    double x0 = 2 * sqrt(theta) / sigma;
    double be = x0 / sqrt(dt) / floor(x0 / sqrt(1.5 * dt));
    double bc = x0 / sqrt(dt) / floor(x0 / sqrt(1.5 * dt) + 1);
    b = (fabs(bc - sqrt(1.5)) < fabs(be - sqrt(1.5))) ? bc : be;
}

void HestonModel::construct_vol_tree()
{
    int num_nodes = 2 * steps - 1;
    int mid = num_nodes / 2;
    vol_tree.resize(steps, vector<shared_ptr<VarianceNode>>(num_nodes, nullptr));
    double x0 = 2 * sqrt(V0) / sigma;
    vol_tree[0][mid] = make_shared<VarianceNode>(x0, V0);

    for (int t = 1; t < steps; t++)
    {

        for (int n = 0; n < 2 * t - 1; n++)
        {
            int parentIdx = mid - t + n;
            auto parent = vol_tree[t - 1][parentIdx];
            if (!parent)
                continue;

            double xt = parent->xt;
            int J = floor((1 / xt) * (0.5 * kappa * (4 * theta / (sigma * sigma) - xt * xt) - 0.5) * sqrt(dt) / b + 1 / (b * b));

            double xu = xt + b * (J + 1) * sqrt(dt);
            double xm = xt + b * J * sqrt(dt);
            double xd = xt + b * (J - 1) * sqrt(dt);

            double vu = (xu < threshold) ? (xu * xu * sigma * sigma) / 4 : 0;
            double vm = (xu < threshold) ? (xm * xm * sigma * sigma) / 4 : 0;
            double vd = (xu < threshold) ? (xd * xd * sigma * sigma) / 4 : 0;

            // Assign nodes
            vol_tree[t][parentIdx - 1] = make_shared<VarianceNode>(vu, xu);
            vol_tree[t][parentIdx] = make_shared<VarianceNode>(vm, xm);
            vol_tree[t][parentIdx + 1] = make_shared<VarianceNode>(vd, xd);

            parent->probs = compute_vol_prob(xt);
        }
    }
}

vector<double> HestonModel::compute_vol_prob(double xt)
{
    vector<double> probs(3, 0.0);
    double muX = (1 / xt) * (0.5 * kappa * (4 * theta / (sigma * sigma) - xt * xt) - 0.5);
    int J = floor(muX * sqrt(dt) / b + 1 / (b * b));

    if (xt > 0)
    {
        probs[0] = 1 / (2 * b * b) - J / 2.0 + 1 / (2 * b) * muX * sqrt(dt); // pu
        probs[1] = 1 - 1 / (b * b);                                          // pm
        probs[2] = 1 / (2 * b * b) + J / 2.0 - 1 / (2 * b) * muX * sqrt(dt); // pd
    }
    else
    {
        double Xu = xt + b * (J + 1) * sqrt(dt);
        double Vu = (Xu * Xu * sigma * sigma) / 4;
        probs[0] = kappa * theta * dt / Vu; // pu
        probs[1] = 0;                       // pm
        probs[2] = 1 - probs[0];            // pd
    }
    return probs;
}

void HestonModel::construct_stock_tree()
{
    int numNodes = 2 * steps - 1;
    stock_tree.resize(steps, vector<shared_ptr<StockNode>>(numNodes, nullptr));

    // Initialize first node
    double h0 = (r - (rho * kappa * theta) / sigma) * dt;
    double Y0 = log(S0) - (rho / sigma) * V0 - h0;

    int mid = numNodes / 2;
    stock_tree[0][mid] = make_shared<StockNode>(Y0, S0);

    // Construct tree based on variance tree
    for (int t = 1; t < steps; ++t)
    {
        int numNodesAtT = 2 * t + 1;
        int startIdx = mid - t;

        for (int n = 0; n < numNodesAtT; ++n)
        {
            int parentIdx = startIdx + n;
            auto parent = stock_tree[t - 1][parentIdx];
            if (!parent)
                continue;

            double Yt = parent->yt;

            // Get corresponding variance value
            auto varianceNode = vol_tree[t - 1][parentIdx];
            if (!varianceNode)
                continue;
            double Vt = varianceNode->vt;

            // Compute node span
            int k_t = 1;
            double sigmaY0 = sqrt(1 - rho * rho) * sqrt(V0);

            // Compute up, middle, down moves
            int I = round(((rho * kappa / sigma - 0.5) * Vt) / (k_t * sigmaY0 * sqrt(dt)));

            double Yu = Yt + (I + 1) * k_t * sigmaY0 * sqrt(dt);
            double Ym = Yt + I * k_t * sigmaY0 * sqrt(dt);
            double Yd = Yt + (I - 1) * k_t * sigmaY0 * sqrt(dt);

            double St_u = exp(Yu + (rho / sigma) * Vt + h0);
            double St_m = exp(Ym + (rho / sigma) * Vt + h0);
            double St_d = exp(Yd + (rho / sigma) * Vt + h0);

            // Assign new nodes
            stock_tree[t][parentIdx - 1] = make_shared<StockNode>(Yu, St_u);
            stock_tree[t][parentIdx] = make_shared<StockNode>(Ym, St_m);
            stock_tree[t][parentIdx + 1] = make_shared<StockNode>(Yd, St_d);

            // Compute probabilities
            vector<double> probs = compute_stock_prob(Vt, Yt);
            parent->probs = probs;
        }
    }
}

vector<double> HestonModel::compute_stock_prob(double Vt, double Yt)
{
    vector<double> probs(3, 0.0);
    int k_t = 1;

    double muY = (rho * kappa / sigma - 0.5) * Vt;
    double sigmaYt = sqrt(1 - rho * rho) * sqrt(Vt);
    double sigmaY0 = sqrt(1 - rho * rho) * sqrt(V0);

    int I = round(muY / (k_t * sigmaY0 * sqrt(dt)));

    double Yu = Yt + (I + 1) * k_t * sigmaY0 * sqrt(dt);
    double Ym = Yt + I * k_t * sigmaY0 * sqrt(dt);
    double Yd = Yt + (I - 1) * k_t * sigmaY0 * sqrt(dt);

    double eu = Yu - Yt - muY * dt;
    double em = Ym - Yt - muY * dt;
    double ed = Yd - Yt - muY * dt;

    probs[0] = 0.5 * (sigmaYt * sigmaYt * dt + em * ed) / (k_t * k_t * sigmaY0 * sigmaY0 * dt);
    probs[1] = -(sigmaYt * sigmaYt * dt + eu * ed) / (k_t * k_t * sigmaY0 * sigmaY0 * dt);
    probs[2] = 0.5 * (sigmaYt * sigmaYt * dt + eu * em) / (k_t * k_t * sigmaY0 * sigmaY0 * dt);

    return probs;
}
