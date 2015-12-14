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

// Line by Line File Reader into String Vector

class StreamLine
{
  std::string data;
public:
  operator std::string() const { return data; }
  friend std::istream & operator>> (std::istream &is, StreamLine &l)
  {
    std::getline(is, l.data);
    l.data.erase(std::remove(l.data.begin(), l.data.end(), '\r'), l.data.end());
    return is;
  }
};
typedef FileData<std::string, StreamLine> LBL_FileReader;

struct String_FileReader 
{
  std::string content;
  String_FileReader (std::string const filename) 
  { 
    std::ifstream t(filename.c_str());
    std::stringstream buffer;
    buffer << t.rdbuf();
    content = buffer.str();
  }
};

