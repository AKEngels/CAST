#ifndef H_ATOM
#define H_ATOM

#include<memory>
#include<string>

namespace scon {
  template<typename Type> class mathmatrix;
}

struct Atom {
  unsigned int serial;
  std::string symbol;
private:
  std::shared_ptr<scon::mathmatrix<double>> center;
};

#endif