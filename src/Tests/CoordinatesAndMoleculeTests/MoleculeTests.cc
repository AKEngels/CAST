#include"MoleculeTests.h"

#include"../../CoordinatesAndMolecules/Molecule.h"
#include"../../Scon/scon_mathmatrix.h"

TEST_F(MoleculeCreatorTests, testMoleculeWithConnectivity) {
  auto creator = std::make_unique<MoleculeCreatorCreatorWithConnectivity>();
  testMolecule(std::move(creator));
}

void MoleculeCreatorTests::testMolecule(std::unique_ptr<MoleculeCreatorCreator> moleculeCreator) {
  EXPECT_TRUE(TestingHelper::checkIfTwoMoleculesAreSame(*moleculeCreator->builMoleculeCreator()->buildMolecule(), expectedMolecule));
}



std::unique_ptr<MoleculeCreator> MoleculeCreatorTests::MoleculeCreatorCreatorWithConnectivity::builMoleculeCreator() const
{
  auto creator = std::make_unique<MoleculeCreatorWithConnectivity>();
  creator->setCoordinates(ExpectedValues::TwoEthanolCoordinates());
  creator->setSymbols(ExpectedValues::TwoEthanolSymbols());
  creator->setConnectivity(ExpectedValues::TwoEthanolConnectivity());
  return creator;
}

std::unique_ptr<MoleculeCreator> MoleculeCreatorTests::MoleculeCreatorCreatorWithoutConnectivity::builMoleculeCreator() const
{
  auto creator = std::make_unique<MoleculeCreatorWithoutConnectivity>();
  creator->setCoordinates(ExpectedValues::TwoEthanolCoordinates());
  creator->setSymbols(ExpectedValues::TwoEthanolSymbols());
  return creator;
}

std::unique_ptr<MoleculeCreator> MoleculeCreatorTests::MoleculeCreatorCreatorMissingInput::builMoleculeCreator() const
{
  return std::make_unique<MoleculeCreatorWithConnectivity>();
}


scon::mathmatrix<double> ExpectedValues::TwoEthanolCoordinates()
{
  return scon::mathmatrix<double>{
    {5.417300, -0.262340, -0.019630},
    { 3.921160, -0.006370, -0.138650 },
    { 5.932260, 0.654480, 0.336910 },
    { 5.605910, -1.090180, 0.696640 },
    { 5.827660, -0.543860, -1.012080 },
    { 3.377400, 0.317230, 1.111260 },
    { 3.716400, 0.802990, -0.874860 },
    { 3.426640, -0.923040, -0.516490 },
    { 3.610670, 1.266010, 1.286850 },
    { 3.010560, -4.008670, -0.019280 },
    { 1.792640, -3.249880, 0.496620 },
    { 2.837990, -4.326430, -1.068770 },
    { 3.911390, -3.359910, 0.015160 },
    { 3.188130, -4.908110, 0.607200 },
    { 1.582130, -2.359650, -0.136690 },
    { 0.902360, -3.910500, 0.440250 },
    { 1.979660, -2.866540, 1.830860 },
  { 2.672720, -2.156050, 1.829110 }, };
}

std::vector<std::string> ExpectedValues::TwoEthanolSymbols()
{
  return std::vector<std::string>{"C",
    "C",
    "H",
    "H",
    "H",
    "O",
    "H",
    "H",
    "H",
    "C",
    "C",
    "H",
    "H",
    "H",
    "H",
    "H",
    "O",
    "H", };
}

ConnectedIndices ExpectedValues::TwoEthanolConnectivity()
{
  return std::vector<std::pair<std::size_t, std::size_t>>({
    {0, 4}, {0,1}, {0,2}, {0,3},
    {1,6}, {1,7},{1,5},
    {4,8},
    {9,11},{9,12},{9,10},{9,13},
    {10,14},{10,15},{10,16},
    {16,17},
    });
}

bool TestingHelper::checkIfTwoMoleculesAreSame(Molecule const& lhs, Molecule const& rhs) {
  auto const& lhsGraph = lhs.getBondGraph();
  auto const& rhsGraph = rhs.getBondGraph();
  auto lhsIteratorPair = boost::vertices(lhsGraph);
  auto rhsIteratorPair = boost::vertices(rhsGraph);

  for (; lhsIteratorPair.first != lhsIteratorPair.second && rhsIteratorPair.first != rhsIteratorPair.second; ++lhsIteratorPair.first, ++rhsIteratorPair.first) {
    if (lhsGraph[*lhsIteratorPair.first].symbol != rhsGraph[*rhsIteratorPair.first].symbol) return false;
    auto lhsAdjacents = boost::adjacent_vertices(*lhsIteratorPair.first, lhsGraph);
    auto rhsAdjacents = boost::adjacent_vertices(*rhsIteratorPair.first, rhsGraph);
    for (; lhsAdjacents.first != lhsAdjacents.second && rhsAdjacents.first != rhsAdjacents.second; ++lhsAdjacents.first, ++rhsAdjacents.first) {
      if (*lhsAdjacents.first != *rhsAdjacents.first) return false;
    }
  }
}

MoleculeCreatorTests::MoleculeCreatorTests() : expectedMolecule(){}
