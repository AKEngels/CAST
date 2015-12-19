/**
CAST 3
matop.h
Purpose: 
Specific algorithms for specific problems are placed here.
Calculations are done using scon::mathmatrix

@author Dustin Kaiser
@version 2.0
*/

#pragma once
#include "scon_mathmatrix.h"
#include "histogram.h"

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
   * Tranforms a (regular) mathmatrix-obj of one frame (xyz, atom_nr)
   * into a mathmatrix-obj suitable for creating the single-trajectory-matrix
   * in cartesian coordinates. This is of course a very specific and not
   * a general transformation.
   */
  Matrix_Class transform_3n_nf(Matrix_Class const& input);

  /**
   * Massweightens a single trajectory matrix, only makes sense in
   * cartesian coordinates, needs additional coords::Coordinates
   * of structure parsed to get atomic masses. Boolean controls wether
   * multiplication with "10e-10" (i.e. conversion of Angstrom
   * to meters as distances) is performed.
   */
  void massweight(Matrix_Class& input, coords::Coordinates const&, bool = true, std::vector<unsigned int> atomsThatAreUsed = std::vector<unsigned int>());

  /**
  * Undos the massweighting of a single trajectory matrix, only makes sense in
  * cartesian coordinates, needs additional coords::Coordinates
  * of structure parsed to get atomic masses. Boolean controls wether
  * multiplication with "10e-10" (i.e. conversion of Angstrom
  * to meters as distances) was performed. See also: void massweight(...)
  */
  void undoMassweight(Matrix_Class& input, coords::Coordinates const&, bool = true, std::vector<unsigned int> atomsThatAreUsed = std::vector<unsigned int>());

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


  /////////////////////////////////////
  //                              /////
  //    E X C L U S I V E L Y     /////
  //            P C A             /////
  //                              /////
  /////////////////////////////////////
  namespace pca
  {
    /**
     * Tranforms a (regular) mathmatrix-obj of one frame (internals(?), atom_nr)
     * into a mathmatrix-obj suitable for creating the single-trajectory-matrix
     * in internal coordinates. This is of course a very specific and not
     * a general transformation.
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    Matrix_Class transform_3n_nf_internal_pca(Matrix_Class const& input);

    /**
     * Tranforms a (regular) mathmatrix-obj of one frame into a mathmatrix-obj
     * suitable for creating the single-trajectory-matrix in PCA task
     * in cartesian coordinates.
     * NOTE: This is a different function because the INPUTFILE options
     * differ although formally the operations in this function are equivalent to
     * the ones in "transform_3n_nf()"
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    Matrix_Class transform_3n_nf_trunc_pca(Matrix_Class const& input);

    /**
     * Preparation for PCA, returns eigenvalues and eigenvectors of covariance matrix
     * as well as pca_modes (columns is frames, rows are modes).
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    void prepare_pca(Matrix_Class const& input, Matrix_Class& eigenvalues, Matrix_Class& eigenvectors, Matrix_Class& pca_modes, int rank = 0);

    /**
     * Writes aproximated probability density of each mode.
     * Histogramming is used to genereate PD.
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    void output_probability_density(Matrix_Class& pca_modes);

    /**
     * Formatted output of PCA-Modes to file
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    void output_pca_modes(Matrix_Class& eigenvalues, Matrix_Class& eigenvectors, Matrix_Class& pca_modes, std::string filename = "pca_modes.dat", std::string additionalInformation = "void");

    /**
     * Reads Eigenvectors and PCA-trajectory from correctly formatted file.
     * File is correctly formatted if it was created using the CAST::matop::output_pca_modes
     * function.
     */
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
     * Performs entropy calculation according to Hnizdo et al.
     * see: (DOI: 10.1002/jcc.20589)
     * If rows > 10 this might be really computationally expensive
     * Basically, many "next-neighbor"-searches are performed using a
     * primitive brute-force approach.
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    float_type hnizdo_wrapper(Matrix_Class const& input);

    /**
     * Performs entropy calculation according to Hnizdo et al.
     * Sums marginal (== 1-dimensional) entropy, i.e. calculating
     * values for each "row" independently and then building sum.
     * see: (DOI: 10.1002/jcc.20589)
     * Marginal Entropy is upper limit to actual entropy which is
     * computationally expensive ( "hnizdo_wrapper()" )
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    float_type hnizdo_m_wrapper(Matrix_Class const& input);

    /**
     * Performs entropy calculation according to Hnizdo et al.
     * Sums marginal (== 1-dimensional) entropy, i.e. calculating
     * values for each "row" independently and then building sum.
     * see: (DOI: 10.1002/jcc.20589)
     * Marginal Entropy is upper limit to actual entropy which is
     * computationally expensive ( "hnizdo_wrapper()" )
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    float_type karplus_wrapper(Matrix_Class const& input);

    /**
     * Performs entropy calculation according to Knapp et al. without corrections
     * for anharmonicity or Mutual Information.
     * Quasi-Harmonic-Approximation used.
     * This method is almost identical to Schlitter's approach.
     * Gives a strict upper limit to the actual entropy.
     * see: (Genome Inform. 2007;18:192-205.)
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    float_type knapp_m_wrapper(Matrix_Class const& input);

    /**
     * Performs entropy calculation according to Knapp et al. with corrections
     * for anharmonicity or Mutual Information.
     * Quasi-Harmonic-Approximation used.
     * Hnizdo's entropy estimator is used to calculate the corrections.
     * Gives a strict upper limit to the actual entropy.
     * see: (Genome Inform. 2007;18:192-205.)
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    float_type knapp_wrapper(Matrix_Class const& input);

    /**
     * Performs entropy calculation according to Schlitter
     * Quasi-Harmonic-Approximation used.
     * Gives a strict upper limit to the actual entropy.
     * see: (doi:10.1016/0009-2614(93)89366-P)
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    float_type schlitter_wrapper(Matrix_Class const& input);

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
    float_type knn_distance(Matrix_Class const& input, unsigned int const& dimension_in, unsigned int const& k_in, std::vector<unsigned int>& row_querypts, unsigned int const& col_querypt);

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
    float_type knn_distance(Matrix_Class const& input, unsigned int const& dimension_in, unsigned int const& k_in, unsigned int const& row_querypt, unsigned int const& col_querypt);

    /**
     * Tranforms a (regular) mathmatrix-obj of one frame (internals(?), atom_nr)
     * into a mathmatrix-obj suitable for creating the single-trajectory-matrix
     * in internal coordinates. This is of course a very specific and not
     * a general transformation.
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    Matrix_Class transform_3n_nf_internal_entropy(Matrix_Class const& input);

    /**
     * Tranforms a (regular) mathmatrix-obj of one frame into a mathmatrix-obj
     * suitable for creating the single-trajectory-matrix in Entropy task
     * in cartesian coordinates.
     * NOTE: This is a different function because the INPUTFILE options
     * differ although formally the operations in this function are equivalent to
     * the ones in "transform_3n_nf()"
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    Matrix_Class transform_3n_nf_trunc_entropy(Matrix_Class const& input);
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
     * Returns rotated mathmatrix-obj according to reference mathmatrix-obj "ref"
     * using Kabsch's-Algorithm with singular value decomposition
     */
    Matrix_Class rotated(Matrix_Class const& input, Matrix_Class const& ref);

    /**
     * rotates mathmatrix-obj according to reference mathmatrix-obj "ref"
     * using Kabsch's-Algorithm with singular value decomposition
     */
    void rotate(Matrix_Class& input, Matrix_Class const& ref);

    /**
     * Aligns mathmatrix-obj to geometric Origin with center of geometry
     * DEPRECATED: center of mass is used in current implementaitions (as it should be)
     */
    void align_center_of_geo(Matrix_Class& input);

    /**
     * Returns a 3D-Vector of the Center of Geometry in cartesian coordinates
     * DEPRECATED: center of mass is used in current implementaitions (as it should be)
     * CAUTION WHEN USING!
     */
    std::vector <float_type> center_of_geo(Matrix_Class const& input);

    /**
     * Aligns mathmatrix-obj's center of mass to Geometric Origin
     * ie.: translational alignment
     */
    void align_center_of_mass(Matrix_Class& input, coords::Coordinates const&);

    /**
     * Returns cRMSD-Value of mathmatrix-obj in accordance to "ref" reference mathmatrix-obj
     */
    float_type calc_rmsd(Matrix_Class const& input, Matrix_Class const& ref);

    /**
    * Returns dRMSD-value of the structure.
    */
    float_type drmsd_calc(Matrix_Class const& input, Matrix_Class const& ref);

    /**
     * Calculates Holm&Sanders Distance of the structure
     *
     * ##DIRECTLY READS INPUTFILE! CAUTION WHEN USING!##
     */
    float_type holmsander_calc(Matrix_Class const& input, Matrix_Class const& ref);
  }
}
