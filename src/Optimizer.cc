#include "Optimizer.h"



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
  auto energyDiff = parentOptimizer.currentVariables.systemEnergy - parentOptimizer.oldVariables->systemEnergy;
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

/*class cartesianToInternalNormHelper {
public:
  cartesianToInternalNormHelper();
  coords::float_type getDiffOfCartesianAndInternalNorm(coords::float_type const internalNorm) {

  }
};*/

void Optimizer::optimize(coords::DL_Coordinates<coords::input::formats::pdb> & coords) {

  initializeOptimization(coords);

  for (auto i = 0; i< 100; ++i) {

    evaluateNewCartesianStructure(coords);

    auto cartesianGradients = getInternalGradientsButReturnCartesianOnes(coords);

    applyHessianChange();

    if (ConvergenceCheck{ i + 1,cartesianGradients,*this }()) {
      std::cout << "Converged after " << i + 1 << " steps!\n";
      break;
    }
    //output.to_stream(std::cout);

    setNewToOldVariables();

  }
  /*std::cout << "Final Structure: \n\n";
  output.to_stream(std::cout);
  std::ofstream ofs("Conf_struc.txyz");
  output.to_stream(ofs);*/
}

void Optimizer::initializeOptimization(coords::DL_Coordinates<coords::input::formats::pdb> & coords) {
  setCartesianCoordinatesForGradientCalculation(coords);
  internalCoordinateSystem.guess_hessian(cartesianCoordinates);

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
}

void Optimizer::evaluateNewCartesianStructure(coords::DL_Coordinates<coords::input::formats::pdb> & coords) {
  //TODO own class for handling the internal Step
  //Appropriate Step Finder got to be inserted here
  //auto dq_step = converter.getInternalStep(oldVariables->systemGradients, hessian);

  //
  //std::cout << "U:\n" << del_mat << "\n\n";
  //converter.applyInternalChange(dq_step);

  coords.set_xyz(ic_core::rep3d_bohr_to_ang(cartesianCoordinates));
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
  auto d_gq = currentVariables.systemGradients - oldVariables->systemGradients;
  auto dq = internalCoordinateSystem.calc_diff(cartesianCoordinates, oldVariables->systemCartesianRepresentation);
  auto term1 = (d_gq*d_gq.t()) / (d_gq.t()*dq)(0, 0);
  auto term2 = ((hessian*dq)*(dq.t()*hessian)) / (dq.t()*hessian*dq)(0, 0);
  hessian += term1 - term2;
}

void Optimizer::setNewToOldVariables() {
  oldVariables->systemEnergy = currentVariables.systemEnergy;
  oldVariables->systemCartesianRepresentation = cartesianCoordinates;
  oldVariables->systemGradients = std::move(currentVariables.systemGradients);
}
