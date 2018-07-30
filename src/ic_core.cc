#include "ic_core.h"
#include "ic_util.h"
#include "energy.h"
#include<iterator>

using coords::float_type;

coords::Representation_3D ic_core::rep3d_bohr_to_ang(coords::Representation_3D const& bohr){
  coords::Representation_3D ang;
  ang.reserve(bohr.size());
  for(auto const& b: bohr){
    ang.emplace_back(b * energy::bohr2ang);
  }
  return ang;
}

coords::Representation_3D ic_core::grads_to_bohr(coords::Representation_3D const& grads){
  coords::Representation_3D bohr_grads;
  bohr_grads.reserve(grads.size());
  for(auto const& g: grads){
    bohr_grads.emplace_back(g / energy::Hartree_Bohr2Kcal_MolAng);
  }
  return bohr_grads;
}

std::shared_ptr<InternalCoordinates::Rotator> ic_core::build_rotation(InternalCoordinates::CartesiansForInternalCoordinates & target,
  std::vector<std::size_t> const& index_vec) {
  coords::Representation_3D reference;
  for (auto const & ind : index_vec) {
    reference.emplace_back(target.at(ind));
  }
  return InternalCoordinates::Rotator::buildRotator(target, index_vec);
}

std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> ic_core::system::create_trans_x() const {

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
  for (auto const & indices : subSystemIndices) {
    result.emplace_back(std::make_unique<InternalCoordinates::TranslationX>(indices));
  }
  return result;
}

std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> ic_core::system::create_trans_y() const {

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
  for (auto const & indices : subSystemIndices) {
    result.emplace_back(std::make_unique<InternalCoordinates::TranslationY>(indices));
  }
  return result;
}

std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> ic_core::system::create_trans_z() const {

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
  for (auto const & indices : subSystemIndices) {
    result.emplace_back(std::make_unique<InternalCoordinates::TranslationZ>(indices));
  }
  return result;
}

std::tuple<ic_core::system::InternalVec, ic_core::system::InternalVec, ic_core::system::InternalVec> 
ic_core::system::createRotationABC(std::vector<InternalCoordinates::Rotations> & rotations) {
  ic_core::system::InternalVec resultA, resultB, resultC;
  for (auto & rotation : rotations) {
    resultA.emplace_back(std::move(rotation.rotationA));
    resultB.emplace_back(std::move(rotation.rotationB));
    resultC.emplace_back(std::move(rotation.rotationC));
    registeredRotators.emplace_back(rotation.rotator);
  }
  return { std::move(resultA), std::move(resultB), std::move(resultC) };
}

std::tuple<ic_core::system::InternalVec, ic_core::system::InternalVec, ic_core::system::InternalVec> ic_core::system::create_rotations() {
  std::vector<InternalCoordinates::Rotations> result;
  for(auto const& indices : subSystemIndices){
    result.emplace_back(ic_core::build_rotation(xyz_, indices)->makeRotations());
  }
  return createRotationABC(result);
}

std::vector<std::vector<float_type>> ic_core::system::deriv_vec(){
  std::vector<std::vector<float_type>> result;

  for (auto const& pic : primitive_internals) {
    result.emplace_back(pic->der_vec(xyz_));
  }

  return result;
}

scon::mathmatrix<float_type>& ic_core::system::ic_Bmat(){
  if(!new_B_matrix){
    return B_matrix;
  }
  B_matrix = del_mat.t()*Bmat();
  new_B_matrix = false;
  return B_matrix;
}

scon::mathmatrix<float_type>& ic_core::system::Bmat() {
  if(!new_B_matrix){
    return B_matrix;
  }
  using Mat = scon::mathmatrix<float_type>;

  auto ders = deriv_vec();

  std::size_t n_rows = ders.size(), n_cols = ders.at(0).size();
  B_matrix = Mat(n_rows, n_cols);
  for (std::size_t i{ 0 }; i < n_rows; ++i) {
    B_matrix.set_row(i, Mat::row_from_vec(ders.at(i)));
  }
  new_B_matrix = false;
  return B_matrix;
}

scon::mathmatrix<float_type>& ic_core::system::Gmat(){
  if(!new_G_matrix){
    return G_matrix;
  }
  Bmat();
  G_matrix = B_matrix * B_matrix.t();
  new_G_matrix=false;
  return G_matrix;
}

scon::mathmatrix<float_type>& ic_core::system::ic_Gmat(){
  if(!new_G_matrix){
    return G_matrix;
  }
  ic_Bmat();
  G_matrix = B_matrix * B_matrix.t();
  new_G_matrix=false;
  return G_matrix;
}

scon::mathmatrix<float_type> ic_core::system::guess_hessian() {
  using Mat = scon::mathmatrix<float_type>;
  using scon::cross;
  using scon::dot;
  using scon::len;

  std::vector<float_type> values;
  for (auto const & pic : primitive_internals) {
    values.emplace_back(pic->hessian_guess(xyz_));
  }
  
  return Mat::col_from_vec(values).diagmat();
}

scon::mathmatrix<float_type>& ic_core::system::delocalize_ic_system() {
  using Mat = scon::mathmatrix<float_type>;

  Mat eigval, eigvec;
  std::tie(eigval, eigvec) = Gmat().eigensym(false);

  auto row_index_vec = eigval.sort_idx();
  auto col_index_vec = eigval.find_idx([](float_type const & a) {
    return std::abs(a) > 1e-6;
  });

  del_mat = eigvec.submat(row_index_vec, col_index_vec);
  new_B_matrix = new_G_matrix = true; //B and G got to be calculated for the new ic_system
  return del_mat;
}

scon::mathmatrix<float_type> ic_core::system::calculate_internal_grads(scon::mathmatrix<float_type> const& g) {
  return ic_Gmat().pinv() * ic_Bmat() * g;
}

scon::mathmatrix<float_type>& ic_core::system::initial_delocalized_hessian(){
  hessian = del_mat.t() * guess_hessian() * del_mat;
  return hessian;
}

scon::mathmatrix<float_type> ic_core::system::calc_prims(coords::Representation_3D const& xyz) const{
  std::vector<float_type> primitives;
  primitives.reserve(primitive_internals.size());

  for(auto const & pic : primitive_internals){
    primitives.emplace_back(pic->val(xyz));
  }
  return scon::mathmatrix<float_type>::row_from_vec(primitives);
}


scon::mathmatrix<float_type> atomsNorm(scon::mathmatrix<float_type> const& norm){
  scon::mathmatrix<float_type> mat(norm.rows(),1);
  for(auto i = 0u; i<norm.rows(); ++i){
    mat(i,0) = norm.row(i).norm();
  }
  return mat;
}


std::pair<float_type,float_type> ic_core::system::gradientRmsValAndMax(scon::mathmatrix<float_type> const& grads){
  auto norms = atomsNorm(grads);
  return {norms.rmsd(), norms.max()};
}

std::pair<float_type,float_type> ic_core::system::displacementRmsValAndMaxTwoStructures(coords::Representation_3D const& oldXyz, coords::Representation_3D const& newXyz){
    auto q = ic_rotation::quaternion(oldXyz, newXyz);
  auto U = ic_rotation::form_rot(q.second);

  auto new_xyz_mat = ic_util::Rep3D_to_Mat(newXyz - ic_util::get_mean(newXyz));
  auto old_xyz_mat = ic_util::Rep3D_to_Mat(oldXyz - ic_util::get_mean(oldXyz));
  auto rot = (U*new_xyz_mat.t()).t();

  old_xyz_mat -= rot;
  old_xyz_mat *= -energy::bohr2ang;
  auto norms = atomsNorm(old_xyz_mat);
  
  return {norms.rmsd(), norms.max()};
}

std::pair<float_type,float_type> ic_core::system::displacementRmsValAndMax()const{
  return displacementRmsValAndMaxTwoStructures(oldVariables->systemCartesianRepresentation, xyz_);
}

void ic_core::system::ConvergenceCheck::writeAndCalcEnergyDiffs(){
  auto energyDiff = internalCoordinateSystem.currentVariables.systemEnergy - internalCoordinateSystem.oldVariables->systemEnergy;
  std::cout << "Energy now: " << std::fixed << internalCoordinateSystem.currentVariables.systemEnergy
            << " Energy diff: " << energyDiff <<"\n";
}

void ic_core::system::ConvergenceCheck::writeAndCalcGradientRmsd(){
  cartesianGradients.reshape(-1,3);
  std::tie(gradientRms, gradientMax) = gradientRmsValAndMax(cartesianGradients);
  std::cout << "GRMS Cartesian: " << gradientRms << "\n";
  std::cout << "GRMS Max Val: " << gradientMax << "\n";
}

void ic_core::system::ConvergenceCheck::writeAndCalcDisplacementRmsd(){
  float_type drms, dmax;
  std::tie(displacementRms, displacementMax) = internalCoordinateSystem.displacementRmsValAndMax();
  std::cout << "DRMS Cartesian: " << displacementRms << "\n";
  std::cout << "DRMS Max Val: " << displacementMax << "\n";
}

bool ic_core::system::ConvergenceCheck::checkConvergence() const {
    return energyDiff < threshEnergy && gradientRms < threshGradientRms && gradientMax < threshGradientMax
            && displacementRms < threshDisplacementRms && threshDisplacementMax;
}


bool ic_core::system::ConvergenceCheck::operator()(){
  std::cout << "----------------------------------------------------------\n";
  std::cout << "Step " << step << "\n";
  
  writeAndCalcEnergyDiffs();
  writeAndCalcGradientRmsd();
  writeAndCalcDisplacementRmsd();
  
  std::cout << "----------------------------------------------------------\n";
  
  return checkConvergence();
}

class cartesianToInternalNormHelper{
public:
    cartesianToInternalNormHelper();
    coords::float_type getDiffOfCartesianAndInternalNorm(coords::float_type const internalNorm){
        
    }
};


void ic_core::system::optimize(coords::DL_Coordinates<coords::input::formats::pdb> & coords){
  
  initializeOptimization(coords);
    
  for(auto i = 0; i< 100; ++i){

    evaluateNewCartesianStructure(coords);
    
    auto cartesianGradients = getInternalGradientsButReturnCartesianOnes(coords);

    applyHessianChange();
    
    if(ConvergenceCheck{i+1,cartesianGradients,*this}()){
      std::cout << "Converged after " << i+1 << " steps!\n";
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

void ic_core::system::initializeOptimization(coords::DL_Coordinates<coords::input::formats::pdb> & coords){
  setCartesianCoordinatesForGradientCalculation(coords);
  initial_delocalized_hessian();
  
  prepareOldVariablesPtr(coords);
}

void ic_core::system::setCartesianCoordinatesForGradientCalculation(coords::DL_Coordinates<coords::input::formats::pdb> & coords){
  coords.set_xyz(ic_core::rep3d_bohr_to_ang(xyz_));
}

void ic_core::system::prepareOldVariablesPtr(coords::DL_Coordinates<coords::input::formats::pdb> & coords){
  oldVariables = std::make_unique<SystemVariables>();

  oldVariables->systemEnergy = coords.g()/energy::au2kcal_mol;
  auto cartesianGradients = scon::mathmatrix<coords::float_type>::col_from_vec(ic_util::flatten_c3_vec(
    ic_core::grads_to_bohr(coords.g_xyz())
  ));
  oldVariables->systemCartesianRepresentation = xyz_;
  oldVariables->systemGradients = calculate_internal_grads(cartesianGradients);
}

void ic_core::system::evaluateNewCartesianStructure(coords::DL_Coordinates<coords::input::formats::pdb> & coords){
    auto dq_step = get_internal_step(oldVariables->systemGradients);

    //std::cout << "U:\n" << del_mat << "\n\n";
    apply_internal_change(dq_step);

    coords.set_xyz(ic_core::rep3d_bohr_to_ang(xyz_));
}

scon::mathmatrix<float_type> ic_core::system::getInternalGradientsButReturnCartesianOnes(coords::DL_Coordinates<coords::input::formats::pdb> & coords){
    currentVariables.systemEnergy = coords.g()/energy::au2kcal_mol;
    auto cartesianGradients = scon::mathmatrix<coords::float_type>::col_from_vec(ic_util::flatten_c3_vec(
      ic_core::grads_to_bohr(coords.g_xyz())
    ));

    currentVariables.systemGradients = calculate_internal_grads(cartesianGradients);
    
    return cartesianGradients;
}


void ic_core::system::applyHessianChange(){
    auto d_gq = currentVariables.systemGradients - oldVariables->systemGradients;
    auto dq = calc_diff(xyz_, oldVariables->systemCartesianRepresentation);
    auto term1 = (d_gq*d_gq.t())/(d_gq.t()*dq)(0,0);
    auto term2 = ((hessian*dq)*(dq.t()*hessian))/(dq.t()*hessian*dq)(0,0);
    hessian += term1 - term2;
}

void ic_core::system::setNewToOldVariables(){
    oldVariables->systemEnergy = currentVariables.systemEnergy;
    oldVariables->systemCartesianRepresentation = xyz_;
    oldVariables->systemGradients = std::move(currentVariables.systemGradients);
}