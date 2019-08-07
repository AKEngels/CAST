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
#include <stdexcept>

using float_type = coords::float_type;
using uint_type = std::size_t;

typedef scon::mathmatrix<float_type> Matrix_Class;

namespace matop
{
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