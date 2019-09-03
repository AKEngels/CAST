#include "Optimizer.h"
#include "InternalCoordinates.h"
#include "InternalCoordinateUtilities.h"
#include "BestFitRotation.h"

#include "../Scon/scon_mathmatrix.h"
#include<regex>

struct Optimizer::SystemVariables {
	SystemVariables() : systemEnergy{}, systemGradients{ std::make_unique<scon::mathmatrix<internals::float_type>>() }, internalValues{ std::make_unique<scon::mathmatrix<internals::float_type>>() }, systemCartesianRepresentation{} {}
	internals::float_type systemEnergy;
	std::unique_ptr<scon::mathmatrix<internals::float_type>> systemGradients;
	std::unique_ptr<scon::mathmatrix<internals::float_type>> internalValues;
	CartesianType systemCartesianRepresentation;
};


Optimizer::Optimizer(internals::PrimitiveInternalCoordinates & internals, CartesianType const& cartesians)
	: internalCoordinateSystem{ internals }, cartesianCoordinates{ std::make_unique<CartesianType>(cartesians) },
	converter{ internalCoordinateSystem, *cartesianCoordinates }, hessian{ std::make_unique<scon::mathmatrix<internals::float_type>>(internalCoordinateSystem.guess_hessian(*cartesianCoordinates)) }, trustRadius{ 0.1 }, expectedChangeInEnergy{ 0.0 }, stepSize{ std::make_unique<scon::mathmatrix<internals::float_type>>() }, currentVariables{}, oldVariables{} {}

Optimizer::~Optimizer() = default;

scon::mathmatrix<internals::float_type>& Optimizer::getHessian() { return *hessian; }
scon::mathmatrix<internals::float_type> const& Optimizer::getHessian() const { return *hessian; }
void Optimizer::setHessian(scon::mathmatrix<internals::float_type> && newHessian) { *hessian = std::move(newHessian); }
void Optimizer::setHessian(scon::mathmatrix<internals::float_type> const& newHessian) { *hessian = newHessian; }
InternalCoordinates::CartesiansForInternalCoordinates const& Optimizer::getXyz() const { return *cartesianCoordinates; }



std::pair<internals::float_type, internals::float_type> Optimizer::gradientRmsValAndMax(scon::mathmatrix<internals::float_type> const& grads) {
  //auto norms = atomsNorm(grads);
  return { grads.rmsd(), grads.max() };
}

std::pair<internals::float_type, internals::float_type> Optimizer::displacementRmsValAndMax()const {
  return oldVariables->systemCartesianRepresentation.displacementRmsValAndMaxTwoStructures(*cartesianCoordinates);
}

void Optimizer::ConvergenceCheck::writeAndCalcEnergyDiffs() {
  /*auto */energyDiff = parentOptimizer.currentVariables->systemEnergy - parentOptimizer.oldVariables->systemEnergy;
  std::cout << "Energy now: " << std::fixed << parentOptimizer.currentVariables->systemEnergy
    << " Energy diff: " << energyDiff << "\n";
}

void Optimizer::ConvergenceCheck::writeAndCalcGradientRmsd() {
  //cartesianGradients.reshape(-1, 3);
  std::tie(gradientRms, gradientMax) = gradientRmsValAndMax(projectedGradients);
  std::cout << "GRMS Internal: " << gradientRms << "\n";
  std::cout << "GRMS Max Val: " << gradientMax << "\n";
}

void Optimizer::ConvergenceCheck::writeAndCalcDisplacementRmsd() {
  std::tie(displacementRms, displacementMax) = parentOptimizer.displacementRmsValAndMax();
  std::cout << "DRMS Cartesian: " << displacementRms << "\n";
  std::cout << "DRMS Max Val: " << displacementMax << "\n";
}

bool Optimizer::ConvergenceCheck::checkConvergence() const {
  return energyDiff < threshEnergy && gradientRms < threshGradientRms && gradientMax < threshGradientMax
    && displacementRms < threshDisplacementRms && displacementMax < threshDisplacementMax;
}


bool Optimizer::ConvergenceCheck::operator()() {
  std::cout << "----------------------------------------------------------\n";
  std::cout << "Step " << step << "\n";

  writeAndCalcEnergyDiffs();
  writeAndCalcGradientRmsd();
  writeAndCalcDisplacementRmsd();

  std::cout << "----------------------------------------------------------\n";

  return checkConvergence();
}

namespace detail{
  class Line : public std::string{
    friend std::istream & operator>>(std::istream & stream, Line & line){
      return std::getline(stream, line);
    }
  };
}

std::vector<std::string> getAllContentOfFile(std::string const& fileName){
  std::ifstream ifs(fileName);
  return std::vector<std::string>((std::istream_iterator<detail::Line>(ifs)), std::istream_iterator<detail::Line>());
}

scon::mathmatrix<internals::float_type> readMatrix(std::string const& fileName){
  auto lines = getAllContentOfFile(fileName);
  std::vector<std::vector<double>> values;
  std::regex pattern("[+-]?\\d*\\d*[eE]?[+-]?\\d*");
  auto cols{ 0u };
  auto rows{ lines.size() };
  for(auto & line : lines){
    std::istringstream iss(line);
    std::vector<double> row;
    std::transform(std::sregex_token_iterator(line.begin(), line.end(), pattern), std::sregex_token_iterator(), std::back_inserter(row), [](auto const& matcher){
      return std::stod(matcher);
    });
    if(cols == 0) cols = static_cast<unsigned int>(row.size());
    if(row.size() != cols) throw std::runtime_error("The lines in file " + fileName + " contain a different amount of floating point numbers.\n");
  }
  scon::mathmatrix<internals::float_type> result(rows, cols);
  for(auto i=0u; i < rows; ++i){
    for(auto j=0u; j< cols; ++j){
      result(i,j) = values.at(i).at(j);
    }
  }
  return result;
}

void Optimizer::optimize(coords::Coordinates & coords) {

  initializeOptimization(coords);
  
  coords::output::formats::xyz_cast output(coords);
  std::ofstream initialStream("InitialStructure.xyz");
  output.to_stream(initialStream);

  for (auto i = 0; i< 500; ++i) {
    std::stringstream ss;
    ss << "Struct" << std::setfill('0') << std::setw(5) << i+1u << ".xyz";
    std::ofstream stepStream(ss.str());
    output.to_stream(stepStream);
    
    //std::stringstream strstr;
    //strstr << "CASTHessianStep" << std::setfill('0') << std::setw(5u) << i + 1u << ".dat";
    //std::ofstream hessianOfs(strstr.str());
    //hessianOfs << std::fixed << std::setprecision(15) << hessian;

    //std::stringstream gradientName;
    //gradientName << "CASTGradientStep" << std::setfill('0') << std::setw(5u) << i + 1u << ".dat";
    //std::ofstream gradientOfs(gradientName.str());
    //gradientOfs << std::fixed << std::setprecision(15) << oldVariables->systemGradients;

    evaluateNewCartesianStructure(coords);

    /*auto cartesianGradients =*/ getInternalGradientsButReturnCartesianOnes(coords);
    std::cout << "Trust Radius: " << trustRadius << std::endl;
    if(changeTrustStepIfNeccessary()) {
      std::cout << "Rejected Step" << std::endl;
      resetStep(coords);
      continue;
    }
        
    applyHessianChange();
    
    auto projectedGradient = internalCoordinateSystem.projectorMatrix(*cartesianCoordinates) * (*currentVariables->systemGradients);
    if (ConvergenceCheck{ i + 1, projectedGradient, *this }()) {
      std::cout << "Converged after " << i + 1 << " steps!\n";
      break;
    }
    setNewToOldVariables();
    
  }
  std::ofstream ofs("FinalStruct.xyz");
  output.to_stream(ofs);
}

void Optimizer::initializeOptimization(coords::Coordinates & coords) {
  setCartesianCoordinatesForGradientCalculation(coords);

  prepareOldVariablesPtr(coords);
}

void Optimizer::setCartesianCoordinatesForGradientCalculation(coords::Coordinates & coords) {
  coords.set_xyz(cartesianCoordinates->toAngstrom());
}

void Optimizer::prepareOldVariablesPtr(coords::Coordinates & coords) {
  oldVariables = std::make_unique<SystemVariables>();

  oldVariables->systemEnergy = coords.g() / energy::au2kcal_mol;
  auto cartesianGradients = scon::mathmatrix<internals::float_type>::col_from_vec(ic_util::flatten_c3_vec(
    ic_util::grads_to_bohr(coords.g_xyz())
  ));
  oldVariables->systemCartesianRepresentation = *cartesianCoordinates;
  *oldVariables->systemGradients = converter.calculateInternalGradients(cartesianGradients);
  *oldVariables->internalValues = converter.calculateInternalValues();
}

void Optimizer::resetStep(coords::Coordinates & coords){
  (*cartesianCoordinates) = oldVariables->systemCartesianRepresentation;
  coords.set_xyz(cartesianCoordinates->toAngstrom());
  converter.reset();
}
		  
void Optimizer::evaluateNewCartesianStructure(coords::Coordinates & coords) {
  auto stepFinder = internalCoordinateSystem.constructStepFinder(converter, *oldVariables->systemGradients, *hessian, *cartesianCoordinates);
  
  stepFinder->appropriateStep(trustRadius);
  expectedChangeInEnergy = stepFinder->getSolBestStep();
  *stepSize = stepFinder->getBestStep();

  cartesianCoordinates->setCartesianCoordnates(stepFinder->extractCartesians());
  coords.set_xyz(cartesianCoordinates->toAngstrom());
}

bool Optimizer::changeTrustStepIfNeccessary() {
  auto differenceInEnergy = currentVariables->systemEnergy - oldVariables->systemEnergy;
  auto quality = (differenceInEnergy) / expectedChangeInEnergy;
  std::cout << "Trust: " << std::boolalpha << (trustRadius > thre_rj) <<
	  " Energy: " << std::boolalpha << (currentVariables->systemEnergy > oldVariables->systemEnergy) << 
	  " Last Thingy: " << std::boolalpha << (quality < -10. || true) <<
	  " Quality bad: " << std::boolalpha << (quality < -1.) << "\n";
  std::cout << "E: " << std::setprecision(10) << currentVariables->systemEnergy << 
	  " Eprev: " << std::setprecision(10) << oldVariables->systemEnergy << 
	  " Expect: " << std::setprecision(10) << expectedChangeInEnergy << 
	  " Quality: " << std::setprecision(10) << quality << std::endl;
  
  if(quality > goodQualityThreshold){
    trustRadius = std::min(0.3, trustRadius * std::sqrt(2.));
  }
  else if(quality < -1. && currentVariables->systemEnergy > oldVariables->systemEnergy && trustRadius > thre_rj){
    trustRadius = std::max(0.0012, trustRadius / 2.);
    auto cartesianNorm = displacementRmsValAndMax().first;
    
    trustRadius = std::max(0.0012, std::min(trustRadius, cartesianNorm/2.));
    return true;
  }
  else if (quality < badQualityThreshold) {
    trustRadius = std::max(0.0012, trustRadius / 2.);
  }
  return false;
}

scon::mathmatrix<internals::float_type> Optimizer::getInternalGradientsButReturnCartesianOnes(coords::Coordinates & coords) {
  currentVariables->systemEnergy = coords.g() / energy::au2kcal_mol;
  auto cartesianGradients = scon::mathmatrix<internals::float_type>::col_from_vec(ic_util::flatten_c3_vec(
    ic_util::grads_to_bohr(coords.g_xyz())
  ));

  *currentVariables->systemGradients = converter.calculateInternalGradients(cartesianGradients);

  return cartesianGradients;
}


void Optimizer::applyHessianChange() {
  static auto i = 0u;
  auto d_gq = (*currentVariables->systemGradients) - (*oldVariables->systemGradients);
  auto dq = *stepSize;//internalCoordinateSystem.calc_diff(cartesianCoordinates, oldVariables->systemCartesianRepresentation);
  
  if (Config::get().general.verbosity> 3u)
  {
    std::stringstream hessianSS;
    hessianSS << "CASTHessianStep" << std::setfill('0') << std::setw(5u) << i + 1u << ".dat";
    std::ofstream hessianOfs(hessianSS.str());
    hessianOfs << std::fixed << std::setprecision(15) << *hessian;

    std::stringstream gradientSS;
    gradientSS << "CASTGradientChange" << std::setfill('0') << std::setw(5u) << i +1u << ".dat";
    std::ofstream gradientOfs(gradientSS.str());
    gradientOfs << std::fixed << std::setprecision(15) << d_gq;

    std::stringstream displacementSS;
    displacementSS << "CASTDisplacementChange" << std::setfill('0') << std::setw(5) << i + 1u << ".dat";
    std::ofstream displacementOfs(displacementSS.str());
    displacementOfs << dq;
    
    /*std::stringstream intGradSS;
    intGradSS << "CASTInternalGradient" << std::setfill('0') << std::setw(5) << i + 1u << ".dat";
    std::ofstream intGradOfs(intGradSS.str());
    intGradOfs << currentVariables.systemGradients;
    
    std::stringstream projGradSS;
    projGradSS << "CASTProjectedGradient" << std::setfill('0') << std::setw(5) << i + 1u << ".dat";
    std::ofstream projGradOfs(projGradSS.str());
    projGradOfs << internalCoordinateSystem.projectorMatrix(cartesianCoordinates) * currentVariables.systemGradients;*/
  }
  auto term1 = (d_gq*d_gq.t()) / (d_gq.t()*dq)(0, 0);
  auto term2 = (((*hessian)*dq)*(dq.t()*(*hessian))) / (dq.t()*(*hessian)*dq)(0, 0);
  *hessian += term1 - term2;
  ++i;
}

void Optimizer::setNewToOldVariables() {
  oldVariables->systemEnergy = currentVariables->systemEnergy;
  oldVariables->systemCartesianRepresentation = *cartesianCoordinates;
  oldVariables->systemGradients = std::move(currentVariables->systemGradients);
}

