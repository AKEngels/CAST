#ifdef GOOGLE_MOCK

#ifndef H_MOLECULE_TESTS
#define H_MOLECULE_TESTS

#include<vector>
#include<string>
#include<memory>
#include<gtest/gtest.h>

#include"../../CoordinatesAndMolecules/Molecule.h"

namespace scon {
  template<typename Type> class mathmatrix;
}
class MoleculeCreator;

namespace ExpectedValues {
  scon::mathmatrix<double> TwoEthanolCoordinates();
  std::vector<std::string> TwoEthanolSymbols();
  std::vector<std::pair<std::size_t, std::size_t>> TwoEthanolConnectivity();
}

namespace TestingHelper {
  bool checkIfTwoMoleculesAreSame(Molecule const& lhs, Molecule const& rhs);
}

class TestMolecule : public Molecule {
public:
  TestMolecule() : Molecule(ExpectedValues::TwoEthanolConnectivity(), 18u){
    bondGraph[0].symbol = "C";
    bondGraph[1].symbol = "C";
    bondGraph[2].symbol = "H";
    bondGraph[3].symbol = "H";
    bondGraph[4].symbol = "H";
    bondGraph[5].symbol = "O";
    bondGraph[6].symbol = "H";
    bondGraph[7].symbol = "H";
    bondGraph[8].symbol = "H";
    bondGraph[9].symbol = "C";
    bondGraph[10].symbol = "C";
    bondGraph[11].symbol = "H";
    bondGraph[12].symbol = "H";
    bondGraph[13].symbol = "H";
    bondGraph[14].symbol = "H";
    bondGraph[15].symbol = "H";
    bondGraph[16].symbol = "O";
    bondGraph[17].symbol = "H";
  }
};

class MoleculeCreatorTests : public testing::Test {
public:
  MoleculeCreatorTests();
  class MoleculeCreatorCreator {
  public:
    virtual ~MoleculeCreatorCreator() = default;
    virtual std::unique_ptr<MoleculeCreator> builMoleculeCreator() const = 0;
  };
  class MoleculeCreatorCreatorWithConnectivity : public MoleculeCreatorCreator {
  public:
    std::unique_ptr<MoleculeCreator> builMoleculeCreator() const override;
  };
  class MoleculeCreatorCreatorWithoutConnectivity : public MoleculeCreatorCreator {
  public:
    std::unique_ptr<MoleculeCreator> builMoleculeCreator() const override;
  };
  class MoleculeCreatorCreatorMissingInput : public MoleculeCreatorCreator {
  public:
    std::unique_ptr<MoleculeCreator> builMoleculeCreator() const override;
  };

  void testMolecule(std::unique_ptr<MoleculeCreatorCreator> moleculeCreator);

private:
  TestMolecule expectedMolecule;
};

#endif
#endif