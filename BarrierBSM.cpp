#include "BarrierBSM.h"
#include "util.h"

double BarrierBSM::compute_d_BS(double A, double B, double vol) const
{
    return (std::log(A / B) + (interest_rate_ + 0.5 * vol * vol) * time_to_maturity_) /
           (vol * std::sqrt(time_to_maturity_));
}

double BarrierBSM::black_scholes_price(double vol, OptionType option_type, double K, double S) const
{
    BSM bsm(K, S, time_to_maturity_, interest_rate_, option_type);
    return bsm(vol);
}

double BarrierBSM::barrier_call_price(double vol, BarrierOptionType barrier_type, double barrier_level) const
{

    double H = barrier_level;
    double S = spot_;
    double K = strike_;
    double r = interest_rate_;
    double T = time_to_maturity_;

    double d_HS = compute_d_BS(H, S, vol);
    double d_SH = compute_d_BS(S, H, vol);

    double discount = std::exp(-r * T);

    double exponent = (2 * (r - 0.5 * vol * vol)) / (vol * vol);
    double scaling_factor = (H > S) ? std::pow(H / S, exponent) : std::pow(S / H, -exponent);

    double H2_over_S = (H * H) / S;
    double CBS_S_K = black_scholes_price(vol, OptionType::EuropeanCall, K, S);
    double CBS_S_H = black_scholes_price(vol, OptionType::EuropeanCall, H, S);
    double CSB_H2_over_S_K = black_scholes_price(vol, OptionType::EuropeanCall, K, H2_over_S);
    double CSB_H2_over_S_H = black_scholes_price(vol, OptionType::EuropeanCall, H, H2_over_S);
    double PBS_H2_over_S_K = black_scholes_price(vol, OptionType::EuropeanPut, K, H2_over_S);
    double PBS_H2_over_S_H = black_scholes_price(vol, OptionType::EuropeanPut, H, H2_over_S);

    double price = 0.0;

    switch (barrier_type)
    {
    case BarrierOptionType::UpAndIn:
        price = scaling_factor *
                    (PBS_H2_over_S_K - PBS_H2_over_S_H + (H - K) * discount * norm_cdf(-d_HS)) +
                CBS_S_H + (H - K) * discount * norm_cdf(d_SH);
        break;
    case BarrierOptionType::UpAndOut:
        price = CBS_S_K - CBS_S_H - (H - K) * discount * norm_cdf(d_SH) -
                scaling_factor * (CSB_H2_over_S_K - CSB_H2_over_S_H - (H - K) * discount * norm_cdf(d_HS));
        break;
    case BarrierOptionType::DownAndIn:
        price = scaling_factor * CSB_H2_over_S_K;
        break;
    case BarrierOptionType::DownAndOut:
        price = CBS_S_K - scaling_factor * CSB_H2_over_S_K;
        break;
    }

    return std::max(0.0, price);
}