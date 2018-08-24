/**
CAST 3
Optimizer.h
Purpose: Definition of the Optimizer for internal coordinates


@author Julian Erdmannsdörfer
@version 3.0
*/

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "coords.h"
#include "ic_core.h"
#include "PrimitiveInternalCoordinates.h"
#include "scon_mathmatrix.h"

class Optimizer {
protected:
  using CartesianType = InternalCoordinates::CartesiansForInternalCoordinates;
public:

  Optimizer(internals::PrimitiveInternalCoordinates & internals, CartesianType const& cartesians) 
	  : internalCoordinateSystem{ internals }, cartesianCoordinates{cartesians},
    converter{ internalCoordinateSystem, cartesianCoordinates }, hessian{ internalCoordinateSystem.guess_hessian(cartesianCoordinates) }, trustRadius{ 0.1 }, expectedChangeInEnergy{ 0.0 } {}

  void optimize(coords::DL_Coordinates<coords::input::formats::pdb> & coords);//To Test

  scon::mathmatrix<coords::float_type>& getHessian() { return hessian; }
  scon::mathmatrix<coords::float_type> const& getHessian() const { return hessian; };
  template<typename Hessian>
  typename std::enable_if<std::is_same<Hessian, scon::mathmatrix<coords::float_type>>::value>::type
    setHessian(Hessian && newHessian) { hessian = std::forward<Hessian>(newHessian); }
  InternalCoordinates::CartesiansForInternalCoordinates const& getXyz() const { return cartesianCoordinates; }

  static scon::mathmatrix<coords::float_type> atomsNorm(scon::mathmatrix<coords::float_type> const& norm);
  static std::pair<coords::float_type, coords::float_type> gradientRmsValAndMax(scon::mathmatrix<coords::float_type> const& grads);
  static std::pair<coords::float_type, coords::float_type> displacementRmsValAndMaxTwoStructures(coords::Representation_3D const& oldXyz, coords::Representation_3D const& newXyz);
protected:
  void initializeOptimization(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
  void setCartesianCoordinatesForGradientCalculation(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
  void prepareOldVariablesPtr(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
  void evaluateNewCartesianStructure(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
  void changeTrustStepIfNeccessary();
  void applyHessianChange();
  void setNewToOldVariables();
  scon::mathmatrix<coords::float_type> getInternalGradientsButReturnCartesianOnes(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
  
  
  internals::PrimitiveInternalCoordinates & internalCoordinateSystem;
  CartesianType cartesianCoordinates;
  internals::InternalToCartesianConverter converter;
  scon::mathmatrix<coords::float_type> hessian;
  coords::float_type trustRadius;
  coords::float_type expectedChangeInEnergy;

  static auto constexpr thre_rj = 0.01;
  static auto constexpr badQualityThreshold = 0.25;
  static auto constexpr goodQualityThreshold = 0.75;
  
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

  std::pair<coords::float_type, coords::float_type> displacementRmsValAndMax()const;

  struct SystemVariables {
    coords::float_type systemEnergy;
    scon::mathmatrix<coords::float_type> systemGradients;
    CartesianType systemCartesianRepresentation;
  };

  SystemVariables currentVariables;
  std::unique_ptr<SystemVariables> oldVariables;
};
#endif
