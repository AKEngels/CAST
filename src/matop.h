/**
CAST 3
matop.h
Purpose: 
Specific algorithms for specific problems are placed here.
Calculations are done using mathmatrix

@author Dustin Kaiser
@version 2.0
*/

#pragma once
#include "scon_mathmatrix.h"
#include "histogram.h"
#include "coords_io.h"
#include <deque>
#include "scon_angle.h"
#include <stdexcept>

typedef mathmatrix<float_type> Matrix_Class;

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

  /////////////////////////////////////
  //                              /////
  //    E X C L U S I V E L Y     /////
  //            P C A             /////
  //                              /////
  /////////////////////////////////////

  /**
  * Reads Eigenvectors and PCA-trajectory from correctly formatted file.
  * File is correctly formatted if it was created using the CAST::matop::output_pca_modes
  * function.
  */
  namespace pca
  {
    void readEigenvectorsAndModes(Matrix_Class& eigenvectors, Matrix_Class& trajectory, std::string& additionalInformation, std::string filename = "pca_modes.dat");
  }

  /////////////////////////////////////
  //                              /////
  //    E X C L U S I V E L Y     /////
  //        E N T R O P Y         /////
  //                              /////
  /////////////////////////////////////
  namespace entropy
  {
    /**
     * Outputs the !SQUARED! next-neighbor distance in eucledean space
     * of the k-nearest neighbor to the query-Point. Rows are dimensions
     * of the data points, columns are the actual data points.
     *
     * @param dimension_in Dimensionality of the search (needs to
     * be identical to row_querypts.size().
     * @param k_in The k-th nearest neighbor will be searched.
     * @param row_querypts std::vector of rows (ie dimensions)
     * used in the search (example {1, 4, 5}).
     * @param col_querypt Columns index of the query Points.
     */
    float_type knn_distance(
      Matrix_Class const& input, size_t 
      const& dimension_in, size_t const& k_in, 
      std::vector<size_t>& row_querypts, 
      size_t const& col_querypt, 
      coords::float_type* buffer = nullptr);

    /**
     * Outputs the !SQUARED! next-neighbor distance in eucledean space
     * of the k-nearest neighbor to the query-Point. Rows are dimensions
     * of the data points, columns are the actual data points.
     *
     * @param dimension_in Dimensionality of the search (will simply use
     * subsequent rows starting from "row_query_Pt").
     * @param k_in The k-th nearest neighbor will be searched.
     * @param row_querypt Row index of the query Point.
     * @param col_querypt Columns index of the query Point.
     */
    float_type knn_distance(Matrix_Class const& input, size_t const& dimension_in, size_t const& k_in, size_t const& row_querypt, size_t const& col_querypt, coords::float_type* buffer = nullptr);
  }


  /////////////////////////////////////
  //                              /////
  //    E X C L U S I V E L Y     /////
  //      A L I G N M E N T       /////
  //                              /////
  /////////////////////////////////////
  namespace align
  {
    /**
     * Returns rotated coords-obj according to reference coords-obj "reference"
     * using Kabsch's-Algorithm
     *
     * @centerOfMassAlignment: Should center-of-mass be aligned prior to
     * rotational alignment? default = true
     */
    coords::Coordinates kabschAligned(coords::Coordinates const& input, coords::Coordinates const& reference, bool centerOfMassAlign = true);

    /**
     * rotates coords-obj according to reference coords-obj "ref"
     * using Kabsch's-Algorithm
     *
     * @param centerOfMassAlignment: Should center-of-mass be aligned prior to rotational alignment? default = true
     * (this calls "centerOfMassAlignment" function before actual alignment procedure)
     */
    void kabschAlignment(coords::Coordinates& input, coords::Coordinates const& reference, bool centerOfMassAlign = true);

    /**
     * Aligns mathmatrix-obj's center of mass to origin of coordinate system
     * ie.: translational alignment
     */
    void centerOfMassAlignment(coords::Coordinates & in);

    /**
     * Returns dRMSD-value of the structure in respect to the reference structure.
     *
     * @param input: Input Structure
     * @param ref: Reference Structure
     */
    float_type drmsd_calc(coords::Coordinates const& input, coords::Coordinates const& ref);

    /**
     * Calculates Holm&Sanders Distance of the structures
     *
     * @param input: Input Structure
     * @param ref: Reference Structure
     * @param holAndSanderDistance: Holm & Sander Distance, unique parameter for this distance metric, see original publication
     *
     */
    float_type holmsander_calc(coords::Coordinates const& input, coords::Coordinates const& ref, double holmAndSanderDistance = 20);
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//                    W R A P P E R F U N C T I O N S                 //
////////////////////////////////////////////////////////////////////////

/**
 * These functions are called from main and perform the operations necessary
 * for their respective task. They read from the Config-(INPUTFILE)-Options.
 */

void alignment(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords);

void pca_proc(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords);