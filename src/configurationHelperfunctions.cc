/**
CAST 3
configurationHelperfunctions.h
Purpose: Functions to parse input / output / etc.
Used mainly by configuration.h / .cc

@author Dustin Kaiser, Susanne Sauer
@version 1.0
*/
#include "configurationHelperfunctions.h"
#include "helperfunctions.h"

namespace config
{
  bool file_exists_readable(std::string const& filename)
  {
    std::ofstream teststream(filename.c_str(), std::ios_base::in);
    return (teststream.is_open() && teststream.good());
  }

  bool bool_from_iss(std::istringstream& in, std::string const& option)
  {
    std::string holder;
    in >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE" || holder == "1")
    {
      return true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE" || holder == "0")
    {
      return false;
    }
    else if (option == std::string())
    {
      throw std::runtime_error("Could not determine boolean input value.");
    }
    else
    {
      throw std::runtime_error("Could not determine boolean input value \"" + holder + "\" for option \"" + option + "\".");
    }
  }

  std::vector<std::size_t> sorted_indices_from_cs_string(std::string str)
  {
    // remove all spaces
    str.erase(std::remove_if(str.begin(), str.end(),
      [](char c) -> bool {return std::isspace(c) > 0; }),
      str.end());
    // replace all commas with spaces
    std::replace(str.begin(), str.end(), ',', ' ');
    std::string d;
    std::stringstream iss(str);
    std::vector<std::size_t> re;
    // get each seperated value as single string d
    while (iss >> d)
    {
      auto dash_pos = d.find('-');
      // if element is a range (contains '-')
      if (dash_pos != std::string::npos)
      {
        d[dash_pos] = ' ';
        std::stringstream pss{ d };
        std::size_t first(0), last(0);
        if (pss >> first && pss >> last)
        {
          if (first <= last)
          {
            for (auto i = first; i <= last; ++i)
            {
              re.push_back(i);
            }
          }
          else
          {
            throw std::runtime_error("Invalid range for indices: '" +
              std::to_string(first) + " - " + std::to_string(last) + "'.");
          }
        }

        else
        {
          throw std::runtime_error("Cannot read index range from '" + d + "'.");
        }
      }
      // throw if non-numeric character is found
      else if (d.find_first_not_of("0123456789") != std::string::npos)
      {
        throw std::runtime_error("Non numeric character found in '" + d + "'.");
      }
      // read number from stringstream of d
      else
      {
        std::stringstream pss{ d };
        std::size_t value;
        if (pss >> value)
        {
          re.push_back(value);
        }
        else
        {
          throw std::runtime_error("Cannot read index from '" + d + "'.");
        }
      }
    }
    // sort resulting numbers
    std::sort(re.begin(), re.end());
    // remove duplicates and return rest
    return std::vector<std::size_t>{re.begin(), std::unique(re.begin(), re.end())};
  }

  std::vector<double> doubles_from_string(std::string str)
  {
    std::vector<double> result;
    std::vector<std::string> stringvec = split(str, ',');
    for (auto i : stringvec)
    {
      if (check_if_number(i) == true) result.emplace_back(std::stod(i));
      else throw std::runtime_error(i + " can't be converted to double.");
    }
    return result;
  }

  std::vector<int> ints_from_string(std::string str)
  {
    if (str == "") return std::vector<int>{};   // if empty string return empty vector
    std::vector<int> result;
    std::vector<std::string> stringvec = split(str, ',');
    for (auto i : stringvec)
    {
      if (check_if_integer(i) == true) result.emplace_back(std::stoi(i));
      else throw std::runtime_error(i + " can't be converted to integer.");
    }
    return result;
  }

}