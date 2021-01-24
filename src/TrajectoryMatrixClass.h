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

private:
  /**
   * Generates matrix representation of
   * a MD trajectory
   */
  void generateCoordinateMatrixfromPCAModesFile(std::string const& filepath, \
    std::size_t start_frame_num = 0u, std::size_t offset = 1u, std::vector<std::size_t> trunc_atoms = std::vector<std::size_t>());

  /**
   * Generates matrix representation of
   * a MD trajectory which has been previously transposed and saved as pca modes.
   */
  void generateCoordinateMatrix(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords, \
    std::size_t start_frame_num = 0u, std::size_t offset = 1u, std::vector<std::size_t> trunc_atoms = std::vector<std::size_t>(), \
    std::size_t io_verbosity = 5u, std::vector<std::size_t> internal_dih = std::vector<std::size_t>(), int ref_frame_alignment = -1, bool massweight = true);

  // This matrix is massweighted when cartesians are used
  Matrix_Class coordsMatrix;
  std::vector<std::size_t> subDims;
};

