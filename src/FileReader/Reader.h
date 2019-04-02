#ifndef H_READER
#define H_READER

#include<vector>
#include<string>
#include<memory>
#include<fstream>

namespace scon {
  template<typename Type> class mathmatrix;
}

class Molecule;

class FileReader {
public:
  virtual void readFile(std::string const& fileName) = 0;
  virtual scon::mathmatrix<double> getCoordinateLines() const = 0;
  virtual std::shared_ptr<Molecule> buildMolecule() = 0;//Not const; May move symbols
};

namespace FileReaderUtil {
  std::string readFile(std::string const& fileName) {
    std::ifstream ifs(fileName);
    return std::string(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
  }
}

#endif