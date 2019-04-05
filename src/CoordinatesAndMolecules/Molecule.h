#ifndef H_MOLECULE
#define H_MOLECULE

#include<memory>
#include<vector>
#include<string>
#include <boost/graph/adjacency_list.hpp>

#include"Atom.h"

namespace scon {
  template<typename Type> class mathmatrix;
}

using ConnectedIndices = std::vector<std::pair<std::size_t, std::size_t>>;

class Molecule;

class MoleculeCreator {
public:
  MoleculeCreator();
  virtual ~MoleculeCreator();
  virtual void setCoordinates(scon::mathmatrix<double>&&);
  virtual void setSymbols(std::vector<std::string>&&);
  virtual std::shared_ptr<Molecule> buildMolecule();

protected:
  struct impl;
  std::unique_ptr<impl> pImpl;
  virtual bool oneFactorIsMissing()const;
};

class MoleculeCreatorWithConnectivity : public MoleculeCreator {
public:
  void setConnectivity(ConnectedIndices&&);
};

class MoleculeCreatorWithoutConnectivity : public MoleculeCreator {
public:
  virtual void setCoordinates(scon::mathmatrix<double>&&) override;
private:
  void findConnectivity();
};

class Molecule {
  using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Atom>;
  friend class MoleculeCreator;
public:
  Graph const& getBondGraph() const { return bondGraph; };
protected:
  Molecule() = default;
  Molecule(ConnectedIndices const&, std::size_t const numberOfAtoms);
  Graph bondGraph;
private:
};

#endif