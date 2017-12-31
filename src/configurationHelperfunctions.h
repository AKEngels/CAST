/**
CAST 3
configurationHelperfunctions.h
Purpose: Functions to parse input / output / etc.
Used mainly by configuration.h / .cc

@author Dustin Kaiser
@version 1.0
*/

#pragma once
#include <sstream>
#include <iostream>
#include <istream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cctype>

namespace config
{
  /*! Checks if file is readable
  *
  * This function tests if a file in the same
  * folder as the CAST executable is readable.
  *
  * @param filename: Full filename of the file
  */
  bool file_exists_readable(std::string const& filename);

  /**
  * Helper function that clips a number to a range
  * Will return either the value if it is inside
  * the boundries set by LOW & HIGH or the respective boundry.
  *
  * @param value: Number to be truncated if neccessary
  * @param LOW: Lower boundry
  * @param HIGH: Upper boundry
  * @return value if LOW < value < HIGH, otherwise return LOW or HIGH
  */
  template<typename T>
  T clip(const T value, const T LOW, const T HIGH)
  {
    using std::min;
    using std::max;
    return min(HIGH, max(LOW, value));
  }

  /*! Creates boolean from integer
   *
   * Helper function that creates a
   * boolean value from an istringstream,
   * scanning for "TRUE", "true", "True" or 1
   * (and equivalently for false)
   *
   * Throws formated output when @value option is given when
   * boolean value cannot be determined.
   */
  bool bool_from_iss(std::istringstream & in, std::string const& option = std::string());


  template<typename T>
  T enum_from_iss(std::istringstream & in)
  {
    int x;
    in >> x;
    return static_cast<T>(x);
  }

  /**
   * Helper function that matches a string
   * to an enum via a sorted "helper-array"
   * If matching could not be performed,
   * enum with value "-1" will be returned.
   *
   * @param arr: Sring Array in same order as enum, used for matching
   * @param ARR_SIZE: Size of the string array
   * @typename ENUM_T: Type of the enum that should be returned
   * @param S: Input String.
   */
  template<typename ENUM_T, std::size_t ARR_SIZE>
  ENUM_T enum_type_from_string_arr(std::string const & S, std::string const (&arr)[ARR_SIZE])
  {
    for (std::size_t i = 0; i < ARR_SIZE; ++i)
    {
      if (S.find(arr[i]) != S.npos) return static_cast<ENUM_T>(i);
    }
    return static_cast<ENUM_T>(-1);
  }

  /*! Creates vector from range of integer-types
  *
  * This creates an std::vector from range based input (istringstream&)
  * only works for integer types
  * EXAMPLE: 3-6 -> 3,4,5,6 (sorted)
  *
  * @todo: remove in favor of "sorted_indices_from_cs_string"
  */
  template<typename T>
  std::vector<T> configuration_range_int(std::vector<std::string>& holder)
  {
    bool none_check = false;
    std::vector<T> temp;

    for (unsigned int i = 0; i < holder.size(); i++)
    {
      if (holder[i] == "none")
      {
        none_check = true;
      }
      bool i_o_done = false;
      while (!i_o_done)
      {
        for (unsigned int j = 0; j < holder[i].size(); j++)
        {
          if (holder[i][j] == *"-")
          {
            int last_comma = -1, next_comma = -1;
            for (int l = j; l >= 0; l--)
            {
              if (holder[i][l] == *",") { last_comma = l; break; }
            }
            for (unsigned int l = j; l < holder[i].size(); l++)
            {
              if (holder[i][l] == *",") { next_comma = l; break; }
              else if (l == (holder[i].size() - 1))
              {
                next_comma = (int)holder[i].size();
                break;
              }
            }
            std::string temp_str;
            int k = stoi(holder[i].substr(last_comma + 1, j - last_comma - 1)) + 1;
            while (k < stoi(holder[i].substr(j + 1, next_comma - j)))
            {
              temp_str.push_back(*",");
              temp_str.append(std::to_string(k));
              k++;
            }
            temp_str.push_back(*",");
            holder[i].erase(j, 1);
            holder[i].insert(j, temp_str);
            i_o_done = false;
            break;
          }
          if (j == holder[i].size() - 1) { i_o_done = true; }
        }
      }
    }
    for (std::size_t i = 0; i < holder.size(); i++)
    {
      if (holder[i] == "none") {
        none_check = true;
        break;
      }
      int last_comma = -1;
      for (std::size_t j = 0; j < holder[i].size(); j++)
      {
        if (holder[i][j] == *",")
        {
          temp.push_back(stoi(holder[i].substr(last_comma + 1, (j - last_comma - 1))));
          last_comma = (int)j;
        }
        else if (j == holder[i].size() - 1)
        {
          temp.push_back(stoi(holder[i].substr(last_comma + 1, (j - last_comma))));
        }
      }
    }
    //Sorting algorithm
    std::size_t size = temp.size();
    std::size_t left = 0;
    while (left < size)
    {
      size_t min = left;
      for (std::size_t i = (left + 1); i < size; i++)
      {
        if (temp[i] < temp[min])
        {
          min = i;
        }
      }

      T holder2 = temp[min];
      temp[min] = temp[left];
      temp[left] = holder2;
      left += 1;
    }
    if (none_check) temp = std::vector<T>{};
    return temp;
  }

  /*! Creates vector from range of float-types
  *
  * This creates an std::vector from range based input (istringstream&)
  * only works for floating types
  * Is much simpler than version for integer types, only packs into vector without filling or sorting
  * EXAMPLE: 3.5,7.8,2 -> 3.5 | 7.8 | 2.0
  */
  template<typename T>
  std::vector<T> configuration_range_float(std::istringstream& cv)
  {
    std::vector<T> temp;
    std::vector<std::string> holder;
    while (cv)
    {
      std::string temp2;
      cv >> temp2;
      holder.push_back(temp2);
    }
    holder.pop_back();
    for (std::size_t i = 0; i < holder.size(); i++)
    {
      temp.push_back(stod(holder[i]));
    }
    return temp;
  }




  /**
   * Helperfunction that (tries to)
   * convert a string to the template
   * type T
   *
   * @param str: Input String
   * @typename T: Class to which the string should be converted
   */
  template<typename T, class CharT, class TraitT, class AllocT>
  T from_string(std::basic_string<CharT, TraitT, AllocT> const &str)
  {
    T tmp;
    std::basic_istringstream<CharT, TraitT, AllocT> is(str);
    is >> tmp;
    return tmp;
  }

  /*! Creates typename from stringstream
   *
   * Helper function that creates a
   * value as specified from an istringstream,
   *
   * Used when reading the config-file.
   */
  template<typename T, class CharT, class TraitT, class AllocT>
  T from_iss(std::basic_istringstream<CharT, TraitT, AllocT> & is)
  {
    T tmp;
    is >> tmp;
    return tmp;
  }

}