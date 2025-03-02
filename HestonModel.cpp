#include "HestonModel.h"

HestonModel::HestonModel(int steps, double S0, double V0, double kappa, double theta,
                         double sigma, double rho, double r, double T, double K, OptionType option_type)
    : steps(steps), S0(S0), V0(V0), kappa(kappa), theta(theta),
      sigma(sigma), rho(rho), r(r), T(T), K(K), option_type(option_type)
{
    dt = T / steps;   // Time step size
    threshold = 1e-8; // Small threshold for variance
    compute_b();      // Compute tree expansion factor `b`
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
    int NR = 2 * steps - 1; // Total number of rows
    int NT = steps;         // Number of time steps
    int M = (NR + 1) / 2;   // Midpoint where volatility starts
    vol_tree.resize(NT, std::vector<std::shared_ptr<VarianceNode>>(NR, nullptr));

    // Initial transformed variance X0
    double X0 = 2 * sqrt(V0) / sigma;

    // Compute `b` parameter (already calculated in compute_b())
    if (b < 1 || b > sqrt(2))
    {
        std::cerr << "Warning: b = " << b << " but should be in range 1 < b < sqrt(2)\n";
    }

    // Initialize the first node (X, V) at time 0
    vol_tree[0][M] = std::make_shared<VarianceNode>(X0, V0);
    vol_tree[0][M]->probs = compute_vol_prob(X0); // Compute probabilities

    // Time 1 node for X (Equations 22 and 30)
    double muX = (1 / X0) * (0.5 * kappa * (4 * theta / (sigma * sigma) - X0 * X0) - 0.5);
    int J = std::floor(muX * sqrt(dt) / b + 1 / (b * b));

    // First nodes at t=1
    for (int i = -1; i <= 1; i++)
    {
        double X_new = X0 + b * (J + i) * sqrt(dt);
        double V_new = (X_new * X_new * sigma * sigma) / 4;
        vol_tree[1][M + i] = std::make_shared<VarianceNode>(X_new, V_new);
        vol_tree[1][M + i]->probs = compute_vol_prob(X_new);
    }

    // Remaining nodes for X
    for (int t = 2; t < NT; t++)
    {
        for (int n = M - t + 1; n < M + t - 1; n++)
        {
            auto parent = vol_tree[t - 1][n];
            if (!parent)
                continue;

            double Xt = parent->xt;
            muX = (1 / Xt) * (0.5 * kappa * (4 * theta / (sigma * sigma) - Xt * Xt) - 0.5);
            J = std::floor(muX * sqrt(dt) / b + 1 / (b * b));

            if (Xt > threshold && Xt * Xt * sigma * sigma / 4 > threshold)
            {
                // Case 1: Nodes where X > 0
                for (int i = -1; i <= 1; i++)
                {
                    double X_new = Xt + b * (J + i) * sqrt(dt);
                    double V_new = (X_new * X_new * sigma * sigma) / 4;
                    vol_tree[t][n + i] = std::make_shared<VarianceNode>(X_new, V_new);
                    vol_tree[t][n + i]->probs = compute_vol_prob(X_new);
                }
            }
            else
            {
                // Case 2: Nodes where X = 0 (Deep Copy)
                if (vol_tree[t - 1][n - 1])
                    vol_tree[t][n - 1] = std::make_shared<VarianceNode>(*vol_tree[t - 1][n - 1]); // Deep copy

                if (vol_tree[t - 1][n])
                    vol_tree[t][n] = std::make_shared<VarianceNode>(*vol_tree[t - 1][n]); // Deep copy

                if (vol_tree[t - 1][n + 1])
                    vol_tree[t][n + 1] = std::make_shared<VarianceNode>(*vol_tree[t - 1][n + 1]); // Deep copy
            }
        }
    }
}

vector<double> HestonModel::compute_vol_prob(double xt)
{
    vector<double> probs(3, 0.0); // [pu, pm, pd]

    // Compute `b` parameter (already calculated in compute_b())
    double X0 = 2 * sqrt(V0) / sigma; // Initial transformed variance

    // Compute `b` using the method from Beliaeva & Nawalkha
    double be = X0 / sqrt(dt) / std::floor(X0 / sqrt(1.5 * dt));
    double bc = X0 / sqrt(dt) / std::floor(X0 / sqrt(1.5 * dt) + 1);
    double b = (fabs(bc - sqrt(1.5)) < fabs(be - sqrt(1.5))) ? bc : be;

    // Compute drift term muX (Equation 22)
    double muX = (1 / xt) * (0.5 * kappa * (4 * theta / (sigma * sigma) - xt * xt) - 0.5);
    int J = std::floor(muX * sqrt(dt) / b + 1 / (b * b));

    if (xt > threshold)
    {
        // Case 1: Probabilities when X > 0 (Equation 28)
        probs[0] = 1 / (2 * b * b) - J / 2.0 + (1 / (2 * b)) * muX * sqrt(dt); // pu
        probs[1] = 1 - 1 / (b * b);                                            // pm
        probs[2] = 1 / (2 * b * b) + J / 2.0 - (1 / (2 * b)) * muX * sqrt(dt); // pd
    }
    else
    {
        // Case 2: Probabilities when X = 0 (Equation 33)
        double Xu = xt + b * (J + 1) * sqrt(dt);
        double Vu = (Xu * Xu * sigma * sigma) / 4;

        probs[0] = (Vu > 0) ? (kappa * theta * dt / Vu) : 1.0; // pu (avoid division by zero)
        probs[1] = 0.0;                                        // pm
        probs[2] = 1 - probs[0];                               // pd
    }

    return probs;
}

void HestonModel::construct_stock_tree()
{
    int numNodes = 2 * steps - 1;
    stock_tree.resize(steps, vector<shared_ptr<StockNode>>(numNodes, nullptr));

    // Initialize the first node at time t=0
    double h0 = (r - (rho * kappa * theta) / sigma) * dt;
    double Y0 = log(S0) - (rho / sigma) * V0 - h0;

    int mid = numNodes / 2;
    stock_tree[0][mid] = make_shared<StockNode>(Y0, S0);

    // Construct tree based on variance tree
    for (int t = 1; t < steps; ++t)
    {
        int numNodesAtT = 2 * t - 1;
        int startIdx = mid - t;

        for (int n = 0; n < numNodesAtT; ++n)
        {
            int parentIdx = startIdx + n;
            auto parent = stock_tree[t - 1][parentIdx];
            if (!parent)
                continue;

            double Yt = parent->yt;

            // Get corresponding variance value from variance tree
            auto varianceNode = vol_tree[t - 1][parentIdx];
            if (!varianceNode)
                continue;
            double Vt = varianceNode->vt;

            // Compute node span k_t based on variance (Equation 11)
            int k_t = std::ceil(sqrt(Vt / V0));
            double sigmaY0 = sqrt(1 - rho * rho) * sqrt(V0);

            // Compute probabilities before assigning nodes
            vector<double> probs = compute_stock_prob(Vt, Yt);
            parent->probs = probs;

            // Compute up, middle, down moves
            int I = round(((rho * kappa / sigma - 0.5) * Vt) / (k_t * sigmaY0 * sqrt(dt)));

            double Yu = Yt + (I + 1) * k_t * sigmaY0 * sqrt(dt);
            double Ym = Yt + I * k_t * sigmaY0 * sqrt(dt);
            double Yd = Yt + (I - 1) * k_t * sigmaY0 * sqrt(dt);

            double St_u = exp(Yu + (rho / sigma) * (Vt - V0) + h0);
            double St_m = exp(Ym + (rho / sigma) * (Vt - V0) + h0);
            double St_d = exp(Yd + (rho / sigma) * (Vt - V0) + h0);

            // Assign new nodes (deep copy to avoid shared memory issues)
            if (!stock_tree[t][parentIdx - 1])
                stock_tree[t][parentIdx - 1] = make_shared<StockNode>(Yu, St_u);
            if (!stock_tree[t][parentIdx])
                stock_tree[t][parentIdx] = make_shared<StockNode>(Ym, St_m);
            if (!stock_tree[t][parentIdx + 1])
                stock_tree[t][parentIdx + 1] = make_shared<StockNode>(Yd, St_d);
        }
    }
}

vector<double> HestonModel::compute_stock_prob(double Vt, double Yt)
{
    vector<double> probs(3, 0.0); // [pu, pm, pd]

    // Compute node span k_t based on variance
    int k_t = (Vt > 0) ? std::ceil(sqrt(Vt / V0)) : 1;

    // Compute sigma_Yt (volatility of transformed stock price)
    double sigmaYt = sqrt(1 - rho * rho) * sqrt(Vt);
    double sigmaY0 = sqrt(1 - rho * rho) * sqrt(V0);

    // Compute drift term for stock price process (Equation 11)
    double muY = (rho * kappa / sigma - 0.5) * Vt;

    // Compute I (scaling factor for movement)
    int I = round(muY / (k_t * sigmaY0 * sqrt(dt)));

    // Compute up, middle, down movements in Yt
    double Yu = Yt + (I + 1) * k_t * sigmaY0 * sqrt(dt);
    double Ym = Yt + I * k_t * sigmaY0 * sqrt(dt);
    double Yd = Yt + (I - 1) * k_t * sigmaY0 * sqrt(dt);

    // Compute errors in movement (Equation 16-17)
    double eu = Yu - Yt - muY * dt;
    double em = Ym - Yt - muY * dt;
    double ed = Yd - Yt - muY * dt;

    // Compute probabilities (Equation 18)
    probs[0] = 0.5 * (sigmaYt * sigmaYt * dt + em * ed) / (k_t * k_t * sigmaY0 * sigmaY0 * dt); // pu
    probs[1] = -(sigmaYt * sigmaYt * dt + eu * ed) / (k_t * k_t * sigmaY0 * sigmaY0 * dt);      // pm
    probs[2] = 0.5 * (sigmaYt * sigmaYt * dt + eu * em) / (k_t * k_t * sigmaY0 * sigmaY0 * dt); // pd

    return probs;
}

double HestonModel::payoff(double s)
{
    if (option_type == OptionType::AmericanCall || option_type == OptionType::EuropeanCall)
    {
        return std::max(s - K, 0.0);
    }
    else
    {
        return std::max(K - s, 0.0);
    }
}

double HestonModel::price_option()
{
    int numNodes = 2 * steps - 1;
    vector<vector<double>> option_tree(steps, vector<double>(numNodes, 0.0));
    int mid = numNodes / 2;

    // Step 1: Compute the payoff at maturity
    for (int i = 0; i < numNodes; i++)
    {
        if (stock_tree[steps - 1][i]) // Ensure the node exists
        {
            double stock_price = stock_tree[steps - 1][i]->st;
            option_tree[steps - 1][i] = payoff(stock_price); // Compute terminal option value
        }
    }

    // Step 2: Backward induction
    for (int t = steps - 2; t >= 0; t--)
    {
        int numNodesAtT = 2 * t - 1;
        int startIdx = mid - t;

        for (int n = 0; n < numNodesAtT; ++n)
        {
            int parentIdx = startIdx + n;
            auto stock_node = stock_tree[t][parentIdx];
            auto variance_node = vol_tree[t][parentIdx];

            if (!stock_node || !variance_node)
                continue;

            double option_price = 0.0;
            vector<double> stock_probs = stock_node->probs;
            vector<double> variance_probs = variance_node->probs;

            // Compute expectation over the 9 child nodes
            int child_idx = 0;
            for (double pS : stock_probs)
            {
                for (double pV : variance_probs)
                {
                    int childPos = parentIdx + (child_idx % 3 - 1) * 3 + (child_idx / 3 - 1);
                    if (childPos >= 0 && childPos < option_tree[t + 1].size())
                    {
                        double jointProb = pS * pV;
                        option_price += jointProb * option_tree[t + 1][childPos];
                    }
                    child_idx++;
                }
            }

            // Discount to present value
            option_tree[t][parentIdx] = exp(-r * dt) * option_price;
        }
    }

    return option_tree[0][mid]; // Return option price at root node
}

void HestonModel::display_vol_tree()
{
    int numNodes = 2 * steps - 1;
    int mid = numNodes / 2;
    std::cout << "===== Variance Tree =====\n";
    for (int t = 0; t < steps; t++)
    {
        std::cout << "Time Step " << t << ":\n";
        int numNodesAtT = t != 0 ? 2 * t - 1 : 1; // Correct node count

        for (int n = 0; n < numNodesAtT; n++)
        {
            int presentIdx = mid - t + n;
            if (presentIdx >= 0 && presentIdx < vol_tree[t].size() && vol_tree[t][presentIdx]) // Ensure valid node
            {
                auto node = vol_tree[t][presentIdx];
                std::cout << std::fixed << std::setprecision(6)
                          << "V[" << t << "][" << presentIdx << "] = " << node->vt
                          << " (Probs: " << node->probs[0] << ", " << node->probs[1] << ", " << node->probs[2] << ")\n";
            }
        }
        std::cout << "------------------------\n";
    }
}

void HestonModel::display_stock_tree()
{
    int numNodes = 2 * steps - 1;
    int mid = numNodes / 2; // Middle index for t=0

    std::cout << "===== Stock Tree =====\n";
    for (int t = 0; t < steps; t++)
    {
        std::cout << "Time Step " << t << ":\n";
        int numNodesAtT = 2 * t - 1; // Correct node count

        for (int n = 0; n < numNodesAtT; n++)
        {
            int presentIdx = mid - t + n; // Adjust index for current step

            if (presentIdx >= 0 && presentIdx < stock_tree[t].size() && stock_tree[t][presentIdx]) // Ensure valid node
            {
                auto node = stock_tree[t][presentIdx];
                std::cout << std::fixed << std::setprecision(6)
                          << "S[" << t << "][" << presentIdx << "] = " << node->st
                          << " (Probs: " << node->probs[0] << ", " << node->probs[1] << ", " << node->probs[2] << ")\n";
            }
        }
        std::cout << "------------------------\n";
    }
}
