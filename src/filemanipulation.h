/**
CAST 3
filemanipulation.h
Purpose:
Reading data smoothly from files.
Usually one does not need to read this.
Instead, just use LBL_FileReader typedef.
Ignorance is bliss...

Caveat: If a file cannot be read, the content of 
"data" obejct will be empty. However, no error will be thrown.

Comments added by Dustin Kaiser October 2016

@author Unknown phd student
@version 1.0
*/

#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <cstddef>
#include <stddef.h>
#include <iterator>
#include <algorithm>
#include <sstream>

namespace filemanip
{
  /**
   * Reads data from an istream, returns a vector with the data.
   * Is not usually manually called but used in FileData struct.
   *
   * @TYPE DATA_T: Type of the data to be stored (example: std::string)
   * @TYPE STREAM_IT_T: Type of the data in stream
   * Caveat: The two above should typically be identical
   */
  template<typename DATA_T, typename STREAM_IT_T>
  std::vector<DATA_T> data_from_stream(std::istream & stream)
  {
    std::vector<DATA_T> data;
    if (stream)
    {
      std::copy(std::istream_iterator<STREAM_IT_T>(stream),
        std::istream_iterator<STREAM_IT_T>(),
        std::back_inserter(data));
    }
    return data;
  }
}

/**
 * Reads data from a file, stores vector with the data in
 * member "data". If the file cannot be read, the content of
 * data will be empty. However, no error will be thrown.
 *
 * @TYPE DATA_T: Type of the data to be stored (example: std::string)
 * @TYPE STREAM_IT_T: Type of the data in stream
 * Caveat: The two above should typically be identical
 * Caveat: If you want to use std::string as type,
 * use predefined typedef "LBL_FileReader".
*/
template<typename DATA_T, typename STREAM_IT_T = DATA_T>
struct FileData
{
  std::vector<DATA_T> data;
  FileData(std::string const & filename)
  {
    std::ifstream input(filename.c_str(), std::ios_base::in);
    data = filemanip::data_from_stream<DATA_T, STREAM_IT_T>(input);
  }
};

/**
 * Dummy object used by FileData struct when reading strings
 * instead of istream. It removes the carriage return character
 * \r from the stream. Specifically tailored for
 * use with FileData struct. Should not be called manually.
 */
class StreamLine
{
  std::string data;
public:
  operator std::string() const { return data; }
  friend std::istream & operator>> (std::istream &is, StreamLine &l)
  {
    std::getline(is, l.data);
    // Remove all \r "carriage return" characters and truncate line accordingly...
    l.data.erase(std::remove(l.data.begin(), l.data.end(), '\r'), l.data.end());
    return is;
  }
};


////////////////////////////
// IMPORTANT TYPEDEF
// Use this if you want to read strings from a file!
typedef FileData<std::string, StreamLine> LBL_FileReader;
////////////////////////////
