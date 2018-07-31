#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "coords.h"
#include "ic_core.h"
#include "PrimitiveInternalCoordinates.h"
#include "scon_mathmatrix.h"

class Optimizer {
public:

  Optimizer(internals::system & internalSystem) : internalCoordinateSystem{internalSystem}{}

  void optimize(coords::DL_Coordinates<coords::input::formats::pdb> & coords);//To Test
protected:
  using CartesianType = InternalCoordinates::CartesiansForInternalCoordinates;
  internals::system & internalCoordinateSystem;
  void initializeOptimization(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
  void setCartesianCoordinatesForGradientCalculation(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
  void prepareOldVariablesPtr(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
  void evaluateNewCartesianStructure(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
  void applyHessianChange();
  void setNewToOldVariables();
  scon::mathmatrix<coords::float_type> getInternalGradientsButReturnCartesianOnes(coords::DL_Coordinates<coords::input::formats::pdb> & coords);

  class ConvergenceCheck {
  public:
    ConvergenceCheck(int step, scon::mathmatrix<coords::float_type> & gxyz, Optimizer const& parent) :
      step{ step },
      cartesianGradients{ gxyz },
      parentOptimizer{ parent },
      energyDiff{ 0.0 },
      gradientRms{ 0.0 },
      displacementRms{ 0.0 },
      gradientMax{ 0.0 },
      displacementMax{ 0.0 } {}

    void writeAndCalcEnergyDiffs();
    void writeAndCalcGradientRmsd();
    void writeAndCalcDisplacementRmsd();
    bool checkConvergence()const;
    bool operator()();
  private:
    int step;
    scon::mathmatrix<coords::float_type> & cartesianGradients;
    Optimizer const& parentOptimizer;

    static auto constexpr threshEnergy = 1.e-6;
    static auto constexpr threshGradientRms = 0.0003;
    static auto constexpr threshDisplacementRms = 0.00045;
    static auto constexpr threshGradientMax = 0.0012;
    static auto constexpr threshDisplacementMax = 0.0018;

    coords::float_type energyDiff;
    coords::float_type gradientRms;
    coords::float_type displacementRms;
    coords::float_type gradientMax;
    coords::float_type displacementMax;
  };

  static scon::mathmatrix<coords::float_type> atomsNorm(scon::mathmatrix<coords::float_type> const& norm);
  static std::pair<coords::float_type, coords::float_type> gradientRmsValAndMax(scon::mathmatrix<coords::float_type> const& grads);
  std::pair<coords::float_type, coords::float_type> displacementRmsValAndMax()const;

  static std::pair<coords::float_type, coords::float_type> displacementRmsValAndMaxTwoStructures(coords::Representation_3D const& oldXyz, coords::Representation_3D const& newXyz);

  struct SystemVariables {
    coords::float_type systemEnergy;
    scon::mathmatrix<coords::float_type> systemGradients;
    CartesianType systemCartesianRepresentation;
  };

  SystemVariables currentVariables;
  std::unique_ptr<SystemVariables> oldVariables;
};

#endif