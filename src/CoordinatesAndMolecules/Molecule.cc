#include"Molecule.h"

#include"../Scon/scon_mathmatrix.h"

struct MoleculeCreator::impl{
  std::unique_ptr<scon::mathmatrix<double>> cartesians;
  std::unique_ptr<std::vector<std::string>> symbols;
  std::unique_ptr<std::vector<ConnectedIndices>> connectivity;
};

bool MoleculeCreator::oneFactorIsMissing() const
{
  return !(pImpl->cartesians) || !(pImpl->symbols) || !(pImpl->connectivity);
}

void MoleculeCreator::setCoordinates(scon::mathmatrix<double>&& cartesiansForMol){
  pImpl->cartesians = std::make_unique<scon::mathmatrix<double>>(std::move(cartesiansForMol));
}

void MoleculeCreator::setSymbols(std::vector<std::string>&& symbolsForMol){
  pImpl->symbols = std::make_unique<std::vector<std::string>>(std::move(symbolsForMol));
}


void MoleculeCreatorWithConnectivity::setConnectivity(std::vector<ConnectedIndices>&& connectivityForMol){
  pImpl->connectivity = std::make_unique<std::vector<ConnectedIndices>>(std::move(connectivityForMol));
}

std::shared_ptr<Molecule> MoleculeCreator::buildMolecule(){
  if (oneFactorIsMissing()) {
    throw std::runtime_error("Could not build the Moleule. Either the symbols, the coordinates, or the connectivity is missing. Thrown in MoleculeCreatorWithoutConnectivity::buildMolecule().");
  }
  auto molecule = std::make_shared<Molecule>(*pImpl->connectivity, pImpl->symbols->size());
  return molecule;
}

void MoleculeCreatorWithoutConnectivity::setCoordinates(scon::mathmatrix<double>&& cartesiansForMol) {
  MoleculeCreator::setCoordinates(std::move(cartesiansForMol));
  findConnectivity();
}

void MoleculeCreatorWithoutConnectivity::findConnectivity(){
}


Molecule::Molecule(ConnectedIndices const & listOfBonds, std::size_t const numberOfAtoms) : bondGraph(listOfBonds.cbegin(), listOfBonds.cend(), numberOfAtoms) {}