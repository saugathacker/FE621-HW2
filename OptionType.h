#pragma once

enum class OptionType
{
    EuropeanCall = 1,
    EuropeanPut = -1,
    AmericanCall = 2,
    AmericanPut = -2
};

enum class BarrierOptionType
{
    UpAndIn = 1,
    UpAndOut = 2,
    DownAndIn = 3,
    DownAndOut = 4
};