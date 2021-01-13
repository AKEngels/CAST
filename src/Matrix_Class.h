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

typedef scon::mathmatrix<coords::float_type> Matrix_Class;

/**
* Class used for representing a TrajectoryMatrix based
* on MD trajectories
* For sample use see task "ENTROPY"
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
    std::size_t io_verbosity = 5, std::vector<std::size_t> internal_dih = std::vector<std::size_t>(), int ref_frame_alignment = -1);

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
    std::size_t io_verbosity = 5u, std::vector<std::size_t> internal_dih = std::vector<std::size_t>(), int ref_frame_alignment = -1);

  // This matrix is massweighted when cartesians are used
  Matrix_Class coordsMatrix;
};


namespace matop
{
  using float_type = ::coords::float_type;
  /////////////////////////////////////
  //                              /////
  // S P E C I F I C   T A S K S  /////
  //                              /////
  /////////////////////////////////////

  /**
   * Transfer mathmatrix-obj to 3Drep, for use in coordinates obj
   * Only makes sense in conventional mathmatrix-obj corresponding to one frame
   * ie.: (xyz, atom_nr)
   */
  coords::Representation_3D transfer_to_3DRepressentation(Matrix_Class const& input);

  /**
   * Constructs mathmatrix-obj from Coordinates as (xyz, atom_nr) (3 x atoms)
   */
  Matrix_Class transform_coordinates(coords::Coordinates& input);

  /**
   * Massweightens a single trajectory matrix, only makes sense in
   * cartesian coordinates, needs additional coords::Coordinates
   * of structure parsed to get atomic masses. Boolean controls wether
   * multiplication with "10e-10" (i.e. conversion of Angstrom
   * to meters as distances) is performed.
   */
  void massweight(Matrix_Class& input, coords::Coordinates const&, bool = true, std::vector<size_t> atomsThatAreUsed = std::vector<size_t>());

  /**
  * Undos the massweighting of a single trajectory matrix, only makes sense in
  * cartesian coordinates, needs additional coords::Coordinates
  * of structure parsed to get atomic masses. Boolean controls wether
  * multiplication with "10e-10" (i.e. conversion of Angstrom
  * to meters as distances) was performed. See also: void massweight(...)
  */
  void undoMassweight(Matrix_Class& input, coords::Coordinates const&, bool = true, std::vector<size_t> atomsThatAreUsed = std::vector<size_t>());

  /**
   * Converts a coords:Coordinates object to a mathmatrix object
   * having the conventional form for single-frame matrices in
   * cartesian coordinates.
   * mathmatrix(xyz, atom_nr)
   */
  Matrix_Class transfer_to_matr(coords::Coordinates const& in);

  /**
   * Converts a coords:Coordinates object to a mathmatrix object
   * having the conventional form for single-frame matrices in
   * internal coordiantes.
   * mathmatrix(dist/angle/dihedral, atom_nr)
   */
  Matrix_Class transfer_to_matr_internal(coords::Coordinates const& in);

  /**
   * Converts a coords:Coordinates object to a mathmatrix object
   * having the conventional form for single-frame matrices in
   * internal coordiantes.
   * mathmatrix(dist/angle/dihedral, atom_nr)
   */
  coords::Representation_Internal transfer_to_internalRepressentation(Matrix_Class const& input);

  /**
   * Tranforms a coords object
   * into a mathmatrix-obj suitable for creating the single-trajectory-matrix (oneliner)
   * This is of course a very specific and not a general transformation.
   */
  Matrix_Class transformToOneline(coords::Coordinates const& coords, std::vector<size_t> const& includedAtoms, bool internalCoordinates = false);
}