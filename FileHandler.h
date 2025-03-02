#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <memory>
#include "Ticker.h" // Include the header where Ticker and OptionData are defined

class FileHandler
{
private:
    char delimiter;
    bool hasHeaders;

public:
    // Default Constructor
    FileHandler(char delimiter = ',', bool hasHeaders = true);

    // Reads CSV and returns a map of tickers
    std::unordered_map<std::string, std::unique_ptr<Ticker>> read_csv_into_ticker_object(const std::string &fileName);
    std::unordered_map<std::string, std::unique_ptr<Ticker>> read_csv_into_ticker_object2(const std::string &fileName);
};
