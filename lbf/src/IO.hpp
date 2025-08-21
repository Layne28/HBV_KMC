#ifndef IO_HPP
#define IO_HPP

#include "ParamDict.hpp"
#include <string_view>
#include <iostream>
#include <fstream>
#include <experimental/filesystem>

void determine_seed(long unsigned int &seed, std::string seed_file, ParamDict &myParams);
bool ends_with(std::string_view str, std::string_view suffix);
#endif