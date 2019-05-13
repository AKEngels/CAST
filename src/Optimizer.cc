#include "Optimizer.h"

#include<regex>


scon::mathmatrix<coords::float_type> Optimizer::atomsNorm(scon::mathmatrix<coords::float_type> const& norm) {
  scon::mathmatrix<coords::float_type> mat(norm.rows(), 1);
  for (auto i = 0u; i<norm.rows(); ++i) {
    mat(i, 0) = norm.row(i).norm();
  }
  return mat;
}

std::pair<coords::float_type, coords::float_type> Optimizer::gradientRmsValAndMax(scon::mathmatrix<coords::float_type> const& grads) {
  auto norms = atomsNorm(grads);
  return { norms.rmsd(), norms.max() };
}

std::pair<coords::float_type, coords::float_type> Optimizer::displacementRmsValAndMaxTwoStructures(coords::Representation_3D const& oldXyz, coords::Representation_3D const& newXyz) {
  auto q = ic_rotation::quaternion(oldXyz, newXyz);
  auto U = ic_rotation::form_rot(q.second);

  auto new_xyz_mat = ic_util::Rep3D_to_Mat(newXyz - ic_util::get_mean(newXyz));
  auto old_xyz_mat = ic_util::Rep3D_to_Mat(oldXyz - ic_util::get_mean(oldXyz));
  auto rot = (U*new_xyz_mat.t()).t();

  old_xyz_mat -= rot;
  old_xyz_mat *= -energy::bohr2ang;
  auto norms = atomsNorm(old_xyz_mat);

  return { norms.rmsd(), norms.max() };
}

std::pair<coords::float_type, coords::float_type> Optimizer::displacementRmsValAndMax()const {
  return displacementRmsValAndMaxTwoStructures(oldVariables->systemCartesianRepresentation, cartesianCoordinates);
}

void Optimizer::ConvergenceCheck::writeAndCalcEnergyDiffs() {
  /*auto */energyDiff = parentOptimizer.currentVariables.systemEnergy - parentOptimizer.oldVariables->systemEnergy;
  std::cout << "Energy now: " << std::fixed << parentOptimizer.currentVariables.systemEnergy
    << " Energy diff: " << energyDiff << "\n";
}

void Optimizer::ConvergenceCheck::writeAndCalcGradientRmsd() {
  cartesianGradients.reshape(-1, 3);
  std::tie(gradientRms, gradientMax) = gradientRmsValAndMax(cartesianGradients);
  std::cout << "GRMS Cartesian: " << gradientRms << "\n";
  std::cout << "GRMS Max Val: " << gradientMax << "\n";
}

void Optimizer::ConvergenceCheck::writeAndCalcDisplacementRmsd() {
  std::tie(displacementRms, displacementMax) = parentOptimizer.displacementRmsValAndMax();
  std::cout << "DRMS Cartesian: " << displacementRms << "\n";
  std::cout << "DRMS Max Val: " << displacementMax << "\n";
}

bool Optimizer::ConvergenceCheck::checkConvergence() const {
  return energyDiff < threshEnergy && gradientRms < threshGradientRms && gradientMax < threshGradientMax
    && displacementRms < threshDisplacementRms && threshDisplacementMax;
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

scon::mathmatrix<coords::float_type> readMatrix(std::string const& fileName){
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
    if(cols == 0) cols = row.size();
    if(row.size() != cols) throw std::runtime_error("The lines in file " + fileName + " contain a different amount of floating point numbers.\n");
  }
  scon::mathmatrix<coords::float_type> result(rows, cols);
  for(auto i=0u; i < rows; ++i){
    for(auto j=0u; j< cols; ++j){
      result(i,j) = values.at(i).at(j);
    }
  }
  return result;
}

void Optimizer::optimize(coords::DL_Coordinates<coords::input::formats::pdb> & coords) {

  initializeOptimization(coords);
  
  coords::output::formats::xyz output(coords);
  std::ofstream initialStream("InitialStucture.xyz");
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

    auto cartesianGradients = getInternalGradientsButReturnCartesianOnes(coords);
    std::cout << "Trust Radius: " << trustRadius << std::endl;
    if(changeTrustStepIfNeccessary()) {
      std::cout << "Rejected Step" << std::endl;
      resetStep(coords);
      continue;
    }
        
    applyHessianChange();

    if (ConvergenceCheck{ i + 1,cartesianGradients,*this }()) {
      std::cout << "Converged after " << i + 1 << " steps!\n";
      break;
    }
    setNewToOldVariables();
    
  }
  std::ofstream ofs("FinalStruct.xyz");
  output.to_stream(ofs);
}

void Optimizer::initializeOptimization(coords::DL_Coordinates<coords::input::formats::pdb> & coords) {
  setCartesianCoordinatesForGradientCalculation(coords);

  prepareOldVariablesPtr(coords);
}

void Optimizer::setCartesianCoordinatesForGradientCalculation(coords::DL_Coordinates<coords::input::formats::pdb> & coords) {
  coords.set_xyz(ic_core::rep3d_bohr_to_ang(cartesianCoordinates));
}

void Optimizer::prepareOldVariablesPtr(coords::DL_Coordinates<coords::input::formats::pdb> & coords) {
  oldVariables = std::make_unique<SystemVariables>();

  oldVariables->systemEnergy = coords.g() / energy::au2kcal_mol;
  auto cartesianGradients = scon::mathmatrix<coords::float_type>::col_from_vec(ic_util::flatten_c3_vec(
    ic_core::grads_to_bohr(coords.g_xyz())
  ));
  oldVariables->systemCartesianRepresentation = cartesianCoordinates;
  oldVariables->systemGradients = converter.calculateInternalGradients(cartesianGradients);
  oldVariables->internalValues = converter.calculateInternalValues();
}

void Optimizer::resetStep(coords::DL_Coordinates<coords::input::formats::pdb> & coords){
  cartesianCoordinates.setCartesianCoordnates(oldVariables->systemCartesianRepresentation);
  coords.set_xyz(ic_core::rep3d_bohr_to_ang(cartesianCoordinates));
  converter.reset();
}
		  
void Optimizer::evaluateNewCartesianStructure(coords::DL_Coordinates<coords::input::formats::pdb> & coords) {
  auto pmat = internalCoordinateSystem.projectorMatrix(cartesianCoordinates);
  auto imat = scon::mathmatrix<coords::float_type>::identity(pmat.rows(), pmat.cols());
  auto projectedHessian = pmat * hessian * pmat + 1000.0 * (imat - pmat);
  internals::AppropriateStepFinder stepFinder(converter, pmat * oldVariables->systemGradients, projectedHessian);
  
  stepFinder.appropriateStep(trustRadius);
  expectedChangeInEnergy = stepFinder.getSolBestStep();
  stepSize = stepFinder.getBestStep();

  cartesianCoordinates.setCartesianCoordnates(stepFinder.extractCartesians());
  coords.set_xyz(ic_core::rep3d_bohr_to_ang(cartesianCoordinates));
}

bool Optimizer::changeTrustStepIfNeccessary() {
  auto differenceInEnergy = currentVariables.systemEnergy - oldVariables->systemEnergy;
  auto quality = (differenceInEnergy) / expectedChangeInEnergy;
  std::cout << "Trust: " << std::boolalpha << (trustRadius > thre_rj) <<
	  " Energy: " << std::boolalpha << (currentVariables.systemEnergy > oldVariables->systemEnergy) << 
	  " Last Thingy: " << std::boolalpha << (quality < -10. || true) <<
	  " Quality bad: " << std::boolalpha << (quality < -1.) << "\n";
  std::cout << "E: " << std::setprecision(10) << currentVariables.systemEnergy << 
	  " Eprev: " << std::setprecision(10) << oldVariables->systemEnergy << 
	  " Expect: " << std::setprecision(10) << expectedChangeInEnergy << 
	  " Quality: " << std::setprecision(10) << quality << std::endl;
  
  if(quality > goodQualityThreshold){
    trustRadius = std::min(0.3, trustRadius * std::sqrt(2.));
  }
  else if(quality < -1. && currentVariables.systemEnergy > oldVariables->systemEnergy && trustRadius > thre_rj){
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

scon::mathmatrix<coords::float_type> Optimizer::getInternalGradientsButReturnCartesianOnes(coords::DL_Coordinates<coords::input::formats::pdb> & coords) {
  currentVariables.systemEnergy = coords.g() / energy::au2kcal_mol;
  auto cartesianGradients = scon::mathmatrix<coords::float_type>::col_from_vec(ic_util::flatten_c3_vec(
    ic_core::grads_to_bohr(coords.g_xyz())
  ));

  currentVariables.systemGradients = converter.calculateInternalGradients(cartesianGradients);

  return cartesianGradients;
}


void Optimizer::applyHessianChange() {
  static auto i = 0u;
  auto d_gq = currentVariables.systemGradients - oldVariables->systemGradients;
  auto dq = stepSize;//internalCoordinateSystem.calc_diff(cartesianCoordinates, oldVariables->systemCartesianRepresentation);
  
  std::stringstream hessianSS;
  hessianSS << "CASTHessianStep" << std::setfill('0') << std::setw(5u) << i + 1u << ".dat";
  std::ofstream hessianOfs(hessianSS.str());
  hessianOfs << std::fixed << std::setprecision(15) << hessian;

  std::stringstream gradientSS;
  gradientSS << "CASTGradientChange" << std::setfill('0') << std::setw(5u) << i +1u << ".dat";
  std::ofstream gradientOfs(gradientSS.str());
  gradientOfs << std::fixed << std::setprecision(15) << d_gq;

  std::stringstream displacementSS;
  displacementSS << "CASTDisplacementChange" << std::setfill('0') << std::setw(5) << i + 1u << ".dat";
  std::ofstream displacementOfs(displacementSS.str());
  displacementOfs << dq;

  auto term1 = (d_gq*d_gq.t()) / (d_gq.t()*dq)(0, 0);
  auto term2 = ((hessian*dq)*(dq.t()*hessian)) / (dq.t()*hessian*dq)(0, 0);
  hessian += term1 - term2;
  ++i;
}

void Optimizer::setNewToOldVariables() {
  oldVariables->systemEnergy = currentVariables.systemEnergy;
  oldVariables->systemCartesianRepresentation = cartesianCoordinates;
  oldVariables->systemGradients = std::move(currentVariables.systemGradients);
}
