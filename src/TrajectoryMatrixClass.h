/**
CAST 3
matop.h
Purpose:
Matrix operations.
For example: Putting coordinates from coords object into a matrix
Transfering coordinates in matrix back to coords object
and so on...
Calculations are done using scon::mathmatrix

@author Dustin Kaiser
@version 2.0
*/

#pragma once
#include "Scon/scon_mathmatrix.h"
#include "coords_io.h"
#include "Scon/scon_angle.h"
#include "PCA.h"
#include <stdexcept>
#include "constants.h"


using uint_type = std::size_t;
using Matrix_Class = pca::Matrix_Class;

/**
* Class used for representing a TrajectoryMatrix based
* on MD trajectories
* For sample use see task "ENTROPY"
* // Rows are DOFs, Columns are frames!
*/
class TrajectoryMatrixRepresentation
{
  using float_type = coords::float_type;
public:
  /**
  * Constructor, subsequently calls (in this order):
  * -> generateCoordinateMatrix
  */
  TrajectoryMatrixRepresentation(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords, \
    std::size_t start_frame_num = 0u, std::size_t offset = 1u, std::vector<std::size_t> trunc_atoms = std::vector<std::size_t>(),\
    std::size_t io_verbosity = 5, std::vector<std::size_t> internal_dih = std::vector<std::size_t>(), int ref_frame_alignment = -1, bool massweight = true);

  /**
  * Constructor, subsequently calls (in this order):
  * -> generateCoordinateMatrixfromPCAModesFile
  */
  TrajectoryMatrixRepresentation(std::string const& filepath, std::size_t start_frame_num = 0u, std::size_t offset = 1u, std::vector<std::size_t> trunc_atoms = std::vector<std::size_t>());

  /**
   * Karplus entropy,
   * basically never works (by design, its a really old approach)
   * see: DOI 10.1021/ma50003a019
   */
  float_type karplus() const;
  /**
  * Performs entropy calculation according to Schlitter
  * Quasi-Harmonic-Approximation used.
  * Gives a strict upper limit to the actual entropy.
  * see: (doi:10.1016/0009-2614(93)89366-P)
  *
  */
  float_type schlitter(float_type const temperatureInKelvin = 300.0) const;

  bool getAreCoordinatesMassWeighted() const
  {
    return this->containsMassWeightedCoordinates;
  }

  Matrix_Class const& getCoordsMatrix(void) const
  {
    return this->coordsMatrix;
  }

  void setCoordsMatrix(Matrix_Class const& in)
  {
    this->coordsMatrix = in;
  }

  std::vector<std::size_t> const& getSubDims() const
  {
    return this->subDims;
  }

  double pcaTransformDraws(Matrix_Class& eigenvaluesPCA, Matrix_Class& eigenvectorsPCA, Matrix_Class& redMasses_out,\
    Matrix_Class& shiftingConstants,\
    Matrix_Class massVector_in, \
    double temperatureInK = 0.0, bool removeDOF = true) const
  {
    //
    // legacy implementation detail
    Matrix_Class drawMatrix = transposed(this->coordsMatrix);
    //
    //
    if (temperatureInK == 0.0 && Config::get().entropy.entropy_temp != 0.0)
    {
      std::cout << "Assuming a temperature of " << Config::get().entropy.entropy_temp << "K for quasiharmonic treatment.\n";
      temperatureInK = Config::get().entropy.entropy_temp;
    }
    //

    std::cout << "Transforming the input coordinates into their PCA modes.\n";
    std::cout << "This directly yields the marginal Quasi - Harmonic - Approx. according to Knapp et. al. without corrections (Genome Inform. 2007; 18:192 - 205)\n";
    if (!this->containsMassWeightedCoordinates)
    {
      std::cout << "ERROR: Mass-Weighted Coordinates are needed as input for the procedure. Aborting..." << std::endl;
      throw std::logic_error("ERROR: Mass-Weighted Coordinates are needed as input for the procedure. Aborting...");
    }
    std::cout << "Commencing calculation..." << std::endl;
    //std::cout << "DEBUG drawMatrix:\n" << drawMatrix << std::endl;
    Matrix_Class input(drawMatrix);
    Matrix_Class cov_matr = Matrix_Class{ input };
    cov_matr = cov_matr - Matrix_Class(input.rows(), input.rows(), 1.) * cov_matr / static_cast<float_type>(input.rows());
    cov_matr = transposed(cov_matr) * cov_matr;
    cov_matr *= (1.f / static_cast<float_type>(input.rows()));
    Matrix_Class eigenvalues;
    Matrix_Class eigenvectors;
    float_type cov_determ = 0.;
    int cov_rank = cov_matr.rank();
    std::tie(eigenvalues, eigenvectors) = cov_matr.eigensym(true);
    // Checks
    if (massVector_in.cols() == 1u)
    {
      if (massVector_in.rows() != eigenvalues.rows())
      {
        throw std::logic_error("Given Mass Vector has wrong dimensionality, aborting.");
        return -1.0;
      }
    }


    if (removeDOF)
    {
      //Remove Eigenvalues that should be zero if cov_matr is singular
      if ((cov_rank < (int)eigenvalues.rows()) || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90))
      {
        std::cout << "Notice: covariance matrix is singular, attempting to fix by truncation of Eigenvalues.\n";
        std::cout << "Details: rank of covariance matrix is " << cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << ".\n";
        size_t temp = std::max(6, int((cov_matr.rows() - cov_rank)));
        eigenvalues.shed_rows(eigenvalues.rows() - temp, eigenvalues.rows() - 1u);
        eigenvectors.shed_cols(eigenvectors.cols() - temp, eigenvectors.cols() - 1u);

      }
      else
      {
        eigenvectors.shed_cols(0, 5);
        eigenvalues.shed_rows(0, 5);
      }
    }



    eigenvaluesPCA = eigenvalues;
    eigenvectorsPCA = eigenvectors;
    Matrix_Class assocRedMasses = calculateReducedMassOfPCAModes(massVector_in, eigenvalues, eigenvectors, this->subDims);
    if (Config::get().general.verbosity > 4)
    {
      std::cout << "eigenvaluesPCA" << std::endl << eigenvalues << std::endl;
      std::cout << "eigenvectorsPCA" << std::endl << eigenvectorsPCA << std::endl;
      std::cout << "DEBUG: assoc red masses:\n" << assocRedMasses << std::endl;
    }
    
    
    Matrix_Class eigenvectors_t(transposed(eigenvectorsPCA));
    Matrix_Class input2(transposed(drawMatrix)); // Root mass weighted cartesian coords, most likely...
    //
    Matrix_Class pcaModes = Matrix_Class(eigenvectors_t * input2);
    const Matrix_Class covarianceMatrixOfPCAModes = pcaModes.covarianceMatrix();
    if (Config::get().general.verbosity > 4)
    {
      std::cout << "DEBUG covariance matrix of mass weighted pca modes:\n";
      for (std::size_t i = 0; i < covarianceMatrixOfPCAModes.rows(); i++)
        std::cout << covarianceMatrixOfPCAModes(i, i) << "\n";
      //std::cout << covarianceMatrixOfPCAModes << std::endl;
      std::cout << std::endl;
    }
    pcaModes = this->unmassweightPCAModes(assocRedMasses, Matrix_Class(eigenvectors_t * input2));
    const Matrix_Class covarianceMatrixOfUnweightedPCAModes = pcaModes.covarianceMatrix();
    if (Config::get().general.verbosity > 4)
    {
      std::cout << "DEBUG covariance matrix of unweighted pca modes:\n";
      for (std::size_t i = 0; i < covarianceMatrixOfUnweightedPCAModes.rows(); i++)
        std::cout << covarianceMatrixOfUnweightedPCAModes(i, i) << "\n";
      //std::cout << covarianceMatrixOfPCAModes << std::endl;
      std::cout << std::endl;
      //std::cout << "PCA-Vec_t:\n" << eigenvectors_t << std::endl;
      //std::cout << "PCA-Modes (unweighted):\n" << this->pcaModes << std::endl; // Nrows are 3xDOFs, NCloumns are Nframes
    }
    const Matrix_Class covarianceMatrixOfINPUT = input2.covarianceMatrix();
    if (Config::get().general.verbosity > 4)
    {
      std::cout << "DEBUG Covariance Matrix of input DrawMatrix:\n";
      for (std::size_t i = 0; i < covarianceMatrixOfINPUT.rows(); i++)
        std::cout << covarianceMatrixOfINPUT(i, i) << "\n";
    }

    //Calculate PCA Frequencies in quasi-harmonic approximation and Entropy in SHO approximation; provides upper limit of entropy
    Matrix_Class pca_frequencies(eigenvalues.rows(), 1u);
    Matrix_Class red_masses(eigenvalues.rows(), 1u);
    Matrix_Class alpha_i(pca_frequencies.rows(), 1u);
    Matrix_Class quantum_entropy(pca_frequencies.rows(), 1u);
    Matrix_Class classical_entropy(pca_frequencies.rows(), 1u);
    Matrix_Class statistical_entropy(pca_frequencies.rows(), 1u);
    //
    // via https://physics.stackexchange.com/questions/401370/normal-modes-how-to-get-reduced-masses-from-displacement-vectors-atomic-masses
    Matrix_Class constant_C(pca_frequencies.rows(), 1u);
    float_type entropy_qho = 0.;
    float_type entropy_cho = 0.;



    for (std::size_t i = 0; i < eigenvalues.rows(); i++)
    {
      if (this->subDims == std::vector<size_t>() || std::find(this->subDims.begin(), this->subDims.end(), i) != this->subDims.end())
      {
        pca_frequencies(i, 0u) = sqrt(constants::boltzmann_constant_kb_SI_units * temperatureInK / eigenvalues(i, 0u));
        if (massVector_in.cols() == 1u && massVector_in.rows() == eigenvalues.rows())
        {
          if (Config::get().general.verbosity > 4)
          {
            std::cout << "....................\n";
            std::cout << "Debug: Mode " << i << std::endl;
            //std::cout << "Debug: kB SI " << constants::boltzmann_constant_kb_SI_units << std::endl;
            std::cout << "Debug: eigenvalues " << eigenvalues(i, 0u) << std::endl;
            std::cout << "Debug: pca_frequencies SI " << pca_frequencies(i, 0u) << std::endl;
            std::cout << "Debug: pca_frequencies cm-1 " << pca_frequencies(i, 0u) / constants::speed_of_light_cm_per_s << std::endl;
            // Assoc red mass of each mode via https://physics.stackexchange.com/questions/401370/normal-modes-how-to-get-reduced-masses-from-displacement-vectors-atomic-masses
          }
          
          double A___normalizationThisEigenvector = 0.;
          double inv_red_mass = 0.0;
          for (std::size_t j = 0; j < eigenvectors.cols(); j++)
          {
            //Each column is one eigenvector
            const double squaredEigenvecValue = eigenvectors(i, j) * eigenvectors(i, j);
            A___normalizationThisEigenvector += squaredEigenvecValue;
            //std::cout << "Debug: A___normalizationThisEigenvector " << A___normalizationThisEigenvector << std::endl;
            const double currentMass = massVector_in(j, 0u);
            //std::cout << "Debug: currentMass " << currentMass << std::endl;
            inv_red_mass += A___normalizationThisEigenvector / currentMass;
            //std::cout << "Debug: inv_red_mass currently  " << inv_red_mass << std::endl;
          }
          //
          red_masses(i, 0u) = 1.0 / inv_red_mass;
          const double& red_mass = red_masses(i, 0u);
          
          //assocRedMasses(i,0u) = red_mass;
          //std::cout << "Debug: Sanity check red_mass: " << assocRedMasses(i, 0u) << std::endl;
          const double squaredStdDev = covarianceMatrixOfPCAModes(i, i);
          //std::cout << "Debug: squaredStdDev (convoluted with red mass) " << squaredStdDev << std::endl;
          //
          const double stdDev_ofPCAMode_inSIUnits = std::sqrt(squaredStdDev) / std::sqrt(red_mass);
          //std::cout << "SDEBUG: sqrt(" << squaredStdDev << ")/sqrt(" << red_mass << ")= " << stdDev_ofPCAMode_inSIUnits << std::endl;
          //const double x_0 = stdDev_ofPCAMode_inSIUnits * std::sqrt(2);
          const double x_0_SI = stdDev_ofPCAMode_inSIUnits * std::sqrt(2);
          const double Sspatial = constants::joules2cal * constants::N_avogadro * (-1.0 * constants::boltzmann_constant_kb_SI_units * (std::log(2 / constants::pi) - std::log(x_0_SI)));
          //std::cout << "Debug: StdDev in SI units " << stdDev_ofPCAMode_inSIUnits << std::endl;
          const double gaussianSSpatial = std::log(stdDev_ofPCAMode_inSIUnits * std::sqrt(2 * constants::pi * std::exp(1.))); //via https://en.wikipedia.org/wiki/Differential_entropy#:~:text=With%20a%20normal%20distribution%2C%20differential,and%20variance%20is%20the%20Gaussian.
          
          const double C1 = constants::boltzmann_constant_kb_SI_units * temperatureInK / constants::h_bar_SI_units / pca_frequencies(i, 0u) * 2. / constants::pi / x_0_SI;
          //std::cout << "Debug: C1 " << C1 << std::endl;
          const double C2 = std::log(C1) + 1.;
          //std::cout << "Debug: C2 " << C2 << std::endl;
          const double C3 = constants::joules2cal * constants::N_avogadro * constants::boltzmann_constant_kb_SI_units * C2;
          if (Config::get().general.verbosity > 4)
          {
            std::cout << "Debug: red_mass " << red_mass << std::endl;
            std::cout << "Debug: Entropy of gaussian with this StdDev " << gaussianSSpatial << std::endl;
            std::cout << "Debug: Sspatial " << Sspatial << std::endl;
            std::cout << "Debug: Sspatial in raw units " << -1.0 * (std::log(2 / constants::pi) - std::log(x_0_SI)) << std::endl;
            std::cout << "Debug: x_0_SI " << x_0_SI << std::endl;
            std::cout << "Debug: C " << C3 << std::endl;
          }
          //C=k*(ln((k*temp)/h_red/freq*2/pi/x_0) + 1)
          //C_dash = C * avogadro
          constant_C(i, 0u) = C3;
        }
        alpha_i(i, 0u) = constants::h_bar_SI_units / (sqrt(constants::boltzmann_constant_kb_SI_units * temperatureInK) * sqrt(eigenvalues(i, 0u)));
        //const double sanityCheck = constants::h_bar_SI_units * pca_frequencies(i, 0u) / constants::boltzmann_constant_kb_SI_units / temperatureInK;
        //std::cout << "Debug: sanitycheck " << sanityCheck << std::endl;
        //These are in units S/k_B (therefore: not multiplied by k_B)
        quantum_entropy(i, 0u) = ((alpha_i(i, 0u) / (exp(alpha_i(i, 0u)) - 1)) - log(1 - exp(-1 * alpha_i(i, 0u)))) * 1.380648813 * 6.02214129 * 0.239005736;
        statistical_entropy(i, 0u) = -1.0 * constants::N_avogadro * constants::boltzmann_constant_kb_SI_units * constants::joules2cal * (log(alpha_i(i, 0u)) -/*this might be plus or minus?!*/ log(sqrt(2. * constants::pi * 2.71828182845904523536)));
        classical_entropy(i, 0u) = -1.0 * constants::N_avogadro * constants::boltzmann_constant_kb_SI_units * constants::joules2cal * (log(alpha_i(i, 0u)) - 1.); // should this be +1??? // The formula written HERE NOW is correct, there is a sign error in the original pape rof Knapp/numata
        if (Config::get().general.verbosity > 4)
        {
          std::cout << "Debug: alpha_i " << alpha_i(i, 0u) << std::endl;
          std::cout << "Debug: classical_entropy " << classical_entropy(i, 0u) << std::endl;
          std::cout << "Debug: statistical_entropy " << statistical_entropy(i, 0u) << std::endl;
          std::cout << "Debug: quantum_entropy " << quantum_entropy(i, 0u) << std::endl;
        }
        //std::cout << "Debug: classical_entropy in raw units " << classical_entropy(i, 0u) / (constants::N_avogadro * constants::boltzmann_constant_kb_SI_units * constants::joules2cal) << std::endl;
        //
        //
        if (!std::isnan(quantum_entropy(i, 0u)))
          entropy_qho += quantum_entropy(i, 0u);
        if (!std::isnan(classical_entropy(i, 0u)))
          entropy_cho += classical_entropy(i, 0u);
      }
    }
    if (Config::get().general.verbosity >= 3u)
    {
      std::cout << "----------\nPrinting PCA_Mode Frequencies and quantum QHA-Entropies:\n";
      std::cout << std::setw(20) << "Mode #" << std::setw(20) << "assoc. red. mass" << std::setw(20) << "cm^-1" << std::setw(20) << "entropy" << std::setw(20) << "% contrib." << std::setw(20) << "C[SI]" << "\n";
      //
      for (std::size_t i = 0; i < eigenvalues.rows(); i++)
      {
        if (this->subDims == std::vector<size_t>() || std::find(this->subDims.begin(), this->subDims.end(), i) != this->subDims.end())
        {
          std::cout << std::setw(20);
          std::cout << i + 1 << std::setw(20) << assocRedMasses(i, 0u) << std::setw(20) << pca_frequencies(i, 0u) / constants::speed_of_light_cm_per_s << std::setw(20) << quantum_entropy(i, 0u);
          std::cout << std::setw(20) << std::to_string(std::round(quantum_entropy(i, 0u) / entropy_qho * 10000) / 100) + " %";
          std::cout << std::setw(20) << constant_C(i, 0u);
          std::cout << "\n";
        }
      }
      std::cout << "----------\n";
    }
    std::cout << "Entropy in quantum QH-approximation from PCA-Modes: " << entropy_qho << " cal / (mol * K)" << std::endl;
    std::cout << "Entropy in classical QH-approximation from PCA-Modes: " << entropy_cho << " cal / (mol * K)" << std::endl;


    shiftingConstants = constant_C;

    redMasses_out = red_masses;
    return entropy_qho;
  }

private:
  /**
   * Generates matrix representation of
   * a MD trajectory
   */
  void generateCoordinateMatrixfromPCAModesFile(std::string const& filepath, \
    std::size_t start_frame_num = 0u, std::size_t offset = 1u, std::vector<std::size_t> trunc_atoms = std::vector<std::size_t>());


  Matrix_Class unmassweightPCAModes(Matrix_Class const& assocRedMasses, Matrix_Class const& pcaModes) const
  {
    if (assocRedMasses.cols() != 1u || pcaModes.rows() != assocRedMasses.rows())
    {
      throw std::logic_error("Cannot un-massweight PCA Modes: Dimensionality of matrices is wrong. Aborting!");
      return Matrix_Class();
    }
    Matrix_Class unweightedPcaModes(pcaModes);
    for (std::size_t i = 0u; i < assocRedMasses.rows(); ++i)
    {
      for (std::size_t j = 0u; j < pcaModes.cols(); ++j)
      {
        const double toBeDivided = unweightedPcaModes(i, j);
        const double divisor = std::sqrt(assocRedMasses(i, 0u));
        unweightedPcaModes(i, j) = toBeDivided / divisor;
      }
    }
    //std::cout << "Debug unweighted:\n" << unweightedPcaModes << "\n";
    return unweightedPcaModes;
  }

  Matrix_Class calculateReducedMassOfPCAModes(Matrix_Class const& massVector, Matrix_Class const& pca_eigenvalues, Matrix_Class const& pca_eigenvectors, std::vector<size_t> const& subDims) const
  {
    Matrix_Class assocRedMasses(pca_eigenvalues.rows(), 1u);
    for (std::size_t i = 0; i < pca_eigenvalues.rows(); i++)
    {
      if (subDims == std::vector<size_t>() || std::find(subDims.begin(), subDims.end(), i) != subDims.end())
      {
        if (massVector.cols() == 1u && massVector.rows() == pca_eigenvalues.rows())
        {
          //std::cout << "....................\n";
          //std::cout << "Debug: Mode " << i << std::endl;
          //std::cout << "Debug: eigenvalues " << pca_eigenvalues << std::endl;
          // Assoc red mass of each mode via https://physics.stackexchange.com/questions/401370/normal-modes-how-to-get-reduced-masses-from-displacement-vectors-atomic-masses
          double A___normalizationThisEigenvector = 0.;
          double inv_red_mass = 0.0;
          for (std::size_t j = 0; j < pca_eigenvectors.cols(); j++)
          {
            //Each column is one eigenvector
            const double squaredEigenvecValue = pca_eigenvectors(i, j) * pca_eigenvectors(i, j);
            A___normalizationThisEigenvector += squaredEigenvecValue;
            //std::cout << "Debug: A___normalizationThisEigenvector " << A___normalizationThisEigenvector << std::endl;
            const double currentMass = massVector(i, 0u);
            //std::cout << "Debug: currentMass " << currentMass << std::endl;
            inv_red_mass += A___normalizationThisEigenvector / currentMass;
            //std::cout << "Debug: inv_red_mass currently  " << inv_red_mass << std::endl;
          }
          //
          const double red_mass = 1.0 / inv_red_mass;
          //std::cout << "Debug: red_mass " << red_mass << std::endl;
          assocRedMasses(i, 0u) = red_mass;
        }
      }
    }
    return assocRedMasses;
  }

  /**
   * Generates matrix representation of
   * a MD trajectory which has been previously transposed and saved as pca modes.
   */
  void generateCoordinateMatrix(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords, \
    std::size_t start_frame_num = 0u, std::size_t offset = 1u, std::vector<std::size_t> trunc_atoms = std::vector<std::size_t>(), \
    std::size_t io_verbosity = 5u, std::vector<std::size_t> internal_dih = std::vector<std::size_t>(), int ref_frame_alignment = -1, bool massweight = true);

  // This matrix is massweighted when cartesians are used
  Matrix_Class coordsMatrix;

  //
  bool containsMassWeightedCoordinates;
  std::vector<std::size_t> subDims;
};

