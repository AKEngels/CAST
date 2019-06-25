#ifndef H_TINKERFILEREADER
#define H_TINKERFILEREADER

#include<regex>

#include"Reader.h"

class TinkerFileReader : public FileReader{
public:
  void readFile(std::string const& fileName) override;
  //scon::mathmatrix<double> getCoordinateLines() const override;
  //std::shared_ptr<Molecule> buildMolecule() override;
};

#endif
