#include "FileHandler.h"

// Constructor implementation
FileHandler::FileHandler(char delimiter, bool hasHeaders)
    : delimiter(delimiter), hasHeaders(hasHeaders) {}

// Function implementation
std::unordered_map<std::string, std::unique_ptr<Ticker>> FileHandler::read_csv_into_ticker_object(const std::string &fileName)
{
    std::unordered_map<std::string, std::unique_ptr<Ticker>> tickers;
    std::ifstream ifile(fileName, std::ios::in);

    if (!ifile.is_open())
    {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return tickers;
    }

    std::string line_;
    if (hasHeaders && std::getline(ifile, line_))
    {
    } // Skip headers if needed

    while (std::getline(ifile, line_))
    {
        std::stringstream ss(line_);
        std::string temp, ticker, expiration, optionType;
        double timeToMaturity, strike, lastPrice, bid, ask, impliedVolatility, spotPrice, interestRate;
        bool inTheMoney;

        std::getline(ss, temp, delimiter); // Ignore first empty column
        std::getline(ss, ticker, delimiter);
        std::getline(ss, expiration, delimiter);
        std::getline(ss, temp, delimiter);
        timeToMaturity = std::stod(temp);
        std::getline(ss, temp, delimiter);
        strike = std::stod(temp);
        std::getline(ss, optionType, delimiter);
        std::getline(ss, temp, delimiter);
        lastPrice = temp.empty() ? 0.0 : std::stod(temp);
        std::getline(ss, temp, delimiter);
        bid = temp.empty() ? 0.0 : std::stod(temp);
        std::getline(ss, temp, delimiter);
        ask = temp.empty() ? 0.0 : std::stod(temp);
        std::getline(ss, temp, delimiter);
        impliedVolatility = temp.empty() ? 0.0 : std::stod(temp);
        std::getline(ss, temp, delimiter);
        inTheMoney = (temp == "True");
        std::getline(ss, temp, delimiter);
        spotPrice = std::stod(temp);
        std::getline(ss, temp, delimiter);
        interestRate = std::stod(temp) / 100.0;

        if (tickers.find(ticker) == tickers.end())
        {
            tickers[ticker] = std::make_unique<Ticker>(ticker, spotPrice, interestRate);
        }

        auto option = std::make_unique<OptionData>(expiration, timeToMaturity, strike, optionType,
                                                   lastPrice, bid, ask, impliedVolatility, inTheMoney);

        tickers[ticker]->addOptionData(std::move(option));
    }

    ifile.close();
    return tickers;
}

// another fucntion for columns matching data from homework1
std::unordered_map<std::string, std::unique_ptr<Ticker>> FileHandler::read_csv_into_ticker_object2(const std::string &fileName)
{
    std::unordered_map<std::string, std::unique_ptr<Ticker>> tickers;
    std::ifstream ifile(fileName, std::ios::in);

    if (!ifile.is_open())
    {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return tickers;
    }

    std::string line_;
    if (hasHeaders && std::getline(ifile, line_))
    {
    } // Skip headers if needed

    while (std::getline(ifile, line_))
    {
        std::stringstream ss(line_);
        std::string temp, ticker, expiration, optionType;
        double timeToMaturity, strike, lastPrice, bid, ask, volume, openInterest, impliedVolatility, spotPrice, interestRate;
        bool inTheMoney;

        std::getline(ss, temp, ','); // Ignore first empty column
        std::getline(ss, ticker, ',');
        std::getline(ss, expiration, ',');
        std::getline(ss, temp, ',');
        timeToMaturity = std::stod(temp);
        std::getline(ss, temp, ',');
        strike = std::stod(temp);
        std::getline(ss, optionType, ',');
        std::getline(ss, temp, ',');
        lastPrice = temp.empty() ? 0.0 : std::stod(temp);
        std::getline(ss, temp, ',');
        bid = temp.empty() ? 0.0 : std::stod(temp);
        std::getline(ss, temp, ',');
        ask = temp.empty() ? 0.0 : std::stod(temp);
        std::getline(ss, temp, ',');
        volume = temp.empty() ? 0.0 : std::stod(temp);
        std::getline(ss, temp, ',');
        openInterest = temp.empty() ? 0.0 : std::stod(temp);
        std::getline(ss, temp, ',');
        impliedVolatility = temp.empty() ? 0.0 : std::stod(temp);
        std::getline(ss, temp, ',');
        inTheMoney = (temp == "True");
        std::getline(ss, temp, ',');
        spotPrice = std::stod(temp);
        std::getline(ss, temp, ',');
        interestRate = std::stod(temp);
        interestRate = interestRate / 100;

        if (tickers.find(ticker) == tickers.end())
        {
            tickers[ticker] = std::make_unique<Ticker>(ticker, spotPrice, interestRate);
        }

        auto option = std::make_unique<OptionData>(expiration, timeToMaturity, strike, optionType,
                                                   lastPrice, bid, ask, impliedVolatility, inTheMoney);

        tickers[ticker]->addOptionData(std::move(option));
    }

    ifile.close();
    return tickers;
}
