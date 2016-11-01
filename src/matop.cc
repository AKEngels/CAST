// For more information see matop.h

#include "matop.h"
namespace matop
{
  /////////////////////////////////////
  //                                ///
  //  S P E C I F I C   T A S K S   ///
  //                                ///
  /////////////////////////////////////

  Matrix_Class transfer_to_matr(coords::Coordinates const& in)
  {
    Matrix_Class out_mat(in.size(), 3u);
    for (size_t l = 0; l < in.size(); l++)
    {
      coords::cartesian_type tempcoord2;
      tempcoord2 = in.xyz(l);
      out_mat(l, 0) = tempcoord2.x();
      out_mat(l, 1) = tempcoord2.y();
      out_mat(l, 2) = tempcoord2.z();
    }
    return transposed(out_mat);
  }

  Matrix_Class transfer_to_matr_internal(coords::Coordinates const& in)
  {
    const size_t sizer = in.size();
    Matrix_Class out_mat(sizer, 3u);
    for (size_t l = 0; l < in.size(); l++)
    {
      out_mat(l, 0) = in.intern(l).radius();
      out_mat(l, 1) = in.intern(l).inclination().radians();
      out_mat(l, 2) = in.intern(l).azimuth().radians();
    }
    return transposed(out_mat);
  }

  Matrix_Class transform_coordinates(coords::Coordinates& input)
  {
	  Matrix_Class output(3, input.size());
	  for (size_t l = 0; l < input.size(); l++)
	  {
		  (output)((l), 0) = input.xyz(l).x();
		  (output)((l), 1) = input.xyz(l).y();
		  (output)((l), 2) = input.xyz(l).z();
	  }
    return output;
  }

  coords::Representation_3D transfer_to_3DRepressentation(Matrix_Class const& input)
  {
    coords::Representation_3D tempcoord1;

    for (size_t i = 0; i < input.cols(); i++)
    {
      coords::Cartesian_Point tempcoord2(input(0u, i), input(1u, i), input(2u, i));
      tempcoord1.push_back(tempcoord2);
    }
    return tempcoord1;
  }

  coords::Representation_Internal transfer_to_internalRepressentation(Matrix_Class const& input)
  {
    coords::Representation_Internal tempcoord1;

    for (size_t i = 0; i < input.cols(); i++)
    {
      coords::internal_type tempcoord2(input(0u, i), scon::ang<float_type>(input(1u, i)), scon::ang<float_type>(input(2u, i)));
      tempcoord1.push_back(tempcoord2);
    }
    return tempcoord1;
  }

  void massweight(Matrix_Class& input, coords::Coordinates const& coords, bool to_meter, std::vector<size_t> atomsThatAreUsed)
  //coords are reference coord object to get atomic masses from
  //boolean controls whether coords should also be multiplied with 10e-10 to convert angstrom to meters)
  {
    if (atomsThatAreUsed.empty())
    {
      for (size_t i = 0; i < input.rows(); i = i + 3)
      {
        double temp = sqrt(coords.atoms(i / 3u).mass() * 1.6605402 * 10e-27);
        if (to_meter)
        {
          temp *= 10e-10;
        }
        for (size_t j = 0; j < input.cols(); j++)
        {
          for (size_t k = 0; k < 3; k++)
          {
            (input)(i + k, j) *= temp;
          }
        }
      }
    }
    else
    {
      for (size_t i = 0; i < input.rows(); i = i + 3)
      {
        double temp = sqrt(coords.atoms(atomsThatAreUsed[i / 3u]).mass() * 1.6605402 * 10e-27);
        if (to_meter)
        {
          temp *= 10e-10;
        }
        for (size_t j = 0; j < input.cols(); j++)
        {
          for (size_t k = 0; k < 3; k++)
          {
            (input)(i + k, j) *= temp;
          }
        }
      }
    }
  }

  void undoMassweight(Matrix_Class& input, coords::Coordinates const& coords, bool to_meter, std::vector<size_t> atomsThatAreUsed)
  {
    if (atomsThatAreUsed.empty())
    {
      for (size_t i = 0; i < input.rows(); i = i + 3)
      {
        double temp = sqrt(coords.atoms(i / 3u).mass() * 1.6605402 * 10e-27);
        if (to_meter)
        {
          temp *= 10e-10;
        }
        for (size_t j = 0; j < input.cols(); j++)
        {
          for (size_t k = 0; k < 3; k++)
          {
            (input)(i + k, j) /= temp;
          }
        }
      }
    }
    else
    {
      for (size_t i = 0; i < input.rows(); i = i + 3)
      {
        double temp = sqrt(coords.atoms(atomsThatAreUsed[i / 3u]).mass() * 1.6605402 * 10e-27);
        if (to_meter)
        {
          temp *= 10e-10;
        }
        for (size_t j = 0; j < input.cols(); j++)
        {
          for (size_t k = 0; k < 3; k++)
          {
            (input)(i + k, j) /= temp;
          }
        }
      }
    }
  }

  Matrix_Class transformToOneline(coords::Coordinates const& coords, std::vector<size_t> const& includedAtoms, bool internalCoordinates)
  {
    //First, some range checks
    if (/*if not all atoms*/ includedAtoms.size() != 0 && includedAtoms[includedAtoms.size() - 1u] > coords.atoms().size())
    {
      throw std::runtime_error("Truncation number is greater than the total number of atoms.");
    }
    else if (internalCoordinates && includedAtoms[0] < 3u)
    {
      throw std::runtime_error("Dihedral with index < 3 specified.");
    }

    if (internalCoordinates)
    {
      // Matrix Size
      Matrix_Class transformed_matrix(1u, (includedAtoms.size() * 2));

      size_t j = 0;
      size_t quicksearch_dih = 0;
      for (size_t i = 0; i < coords.atoms().size(); i++)
      {
        size_t keeper = 0;
        bool checker_dih = false;
        for (size_t l = quicksearch_dih; l < Config::get().PCA.pca_internal_dih.size(); l++)
        {
          if (Config::get().PCA.pca_internal_dih.size() != 0)
          {
            if (Config::get().PCA.pca_internal_dih[l] == i)
            {
              checker_dih = true;
              quicksearch_dih++;
              break;
            }
          }
        }
        if (checker_dih)
        {
          transformed_matrix(0, j + keeper) = cos(coords.intern(i).azimuth().radians());
          transformed_matrix(0, j + keeper + 1) = sin(coords.intern(i).azimuth().radians());
          keeper += 2;
        }
        j = j + keeper;
      }
      return transformed_matrix;
    }
    else
    {
      if (includedAtoms.size() != 0u)
      {
        Matrix_Class transformed_matrix;
        transformed_matrix = Matrix_Class(1u, (3 * includedAtoms.size()));

        int j = 0;
        size_t quicksearch = 0;
        for (size_t i = 0; i < coords.atoms().size(); i++) //iterate over atoms
        {
          bool checker = false;
          for (size_t l = quicksearch; l < includedAtoms.size(); l++) //iterate over vector of atoms to account for
          {
            if (includedAtoms[l] == i)
            {
              checker = true;
              quicksearch++;
              break;
            }
          }
          if (checker)
          {
            transformed_matrix(0, j + 0u) = coords.xyz(i).x();
            transformed_matrix(0, j + 1u) = coords.xyz(i).y();
            transformed_matrix(0, j + 2u) = coords.xyz(i).z();
            j = j + 3u;
          }
        }
        return transformed_matrix;
      }
      else
      {
          Matrix_Class transformed_matrix(1, (coords.atoms().size() * 3u));
          int j = 0;
          for (size_t i = 0; i < coords.atoms().size(); i++)
          {
            for (size_t k = 0; k < 3; k++)
            {
              transformed_matrix(0, j + 0u) = coords.xyz(i).x();
              transformed_matrix(0, j + 1u) = coords.xyz(i).y();
              transformed_matrix(0, j + 2u) = coords.xyz(i).z();
            }
            j = j + 3u;
          }
          return transformed_matrix;
      }
    }
  }

  /////////////////////////////////////
  //                              /////
  //    E X C L U S I V E L Y     /////
  //        E N T R O P Y         /////
  //                              /////
  /////////////////////////////////////
  namespace entropy
  {

    float_type knn_distance(Matrix_Class const& input, size_t const& dimension_in, size_t const& k_in, size_t const& row_queryPt, size_t const& col_queryPt, coords::float_type* buffer)
    {
      float_type temp_distance = 0.0;
      float_type hold_distance;
      float_type *distanceList = buffer;
      if (buffer == nullptr)
      {
        distanceList = new float_type[k_in];
      }
      distanceList[0] = std::numeric_limits<float_type>::max();
      // This iterates over the "n-th" next neighbors
      // (to get the second next neighbor you have to find the first next neighbor etc. )
      for (size_t i = 0; i < k_in; i++)
      {
        // Get max() is initial value for distance comparison.
        // This garantues that the first calcualted value is smaller than
        // hold_distance.
        hold_distance = std::numeric_limits<float_type>::max();

        // This iterates over all points in the set
        for (size_t j = 0; j < input.cols(); j++)
        {
          // Of course we cannot count the distance of an member to itself
          // since it is =0.0
          if (j == col_queryPt) { continue; }

          // For number of dimensions, add the squared distance of the queryPt
          // to the current point ("j") of each dimensions which equals a
          // squared distance in euclidean space
          temp_distance = 0.0;
          for (size_t l = 0; l < dimension_in; l++)
          {
            temp_distance += pow((input)(row_queryPt + l, j) - (input)(row_queryPt + l, col_queryPt), 2);
          }

          // If we are searching for the actual, "first" nearest neighbor
          // we will compare the distance in "hold_distance", initialized as a huge number,
          // to the calcualted temp_distance ("Is this point nearer to the query point
          // than the previously established nearest point?"). We will, of course,
          // keep the smaller (squared) distance of the two.
          if (i == 0 && temp_distance < hold_distance)
          {
            hold_distance = temp_distance;
          }

          // If we are searching for the "(i + 1)-th" nearest neighbot
          // we will compare "hold_distance" to "temp_distance". If temp_distance is smaller, it
          // will be kept, however, only if it is also larger than the previously established
          // "i-th" nearest neighbor. Otherwise we would always just get the absolute nearest neighbor in
          // a set of points, and not the k-th nearest neighbor.
          else if (i > 0 && temp_distance < hold_distance && temp_distance > distanceList[i - 1])
          {
            hold_distance = temp_distance;
          }
        }
        distanceList[i] = hold_distance;
      }
      float_type keeper = distanceList[k_in - 1u];
      if (buffer == nullptr) delete[] distanceList;
      return keeper;
    }

    float_type knn_distance(Matrix_Class const& input, size_t const& dimension_in, size_t const& k_in, std::vector<size_t>& row_queryPts, size_t const& col_queryPt, coords::float_type* buffer)
    //Returns squared distances in the higher-dimensional NN-query case. Needs vector-form input of query Pts
    //Will throw if input is wrong
    {
      //CHECKS
      if (row_queryPts.size() != dimension_in)
      {
        throw("Error in Matrix_Class-function \"knn_distance\": size of \'row_queryPts\' and \'dimension_in\' does not match.");
      }

      float_type temp_distance;
      float_type hold_distance;
      float_type *distanceList = buffer;
      if (buffer == nullptr)
      {
        distanceList = new float_type[k_in];
      }
      distanceList[0] = std::numeric_limits<double>::max();

      // This iterates over the "n-th" next neighbors
      // (to get the second next neighbor you have to find the first next neighbor etc. )
      for (size_t i = 0; i < k_in; i++)
      {
        // Get max() is initial value for distance comparison.
        // This garantues that the first calcualted value is smaller than
        // hold_distance.
        hold_distance = std::numeric_limits<float_type>::max();

        for (size_t j = 0; j < input.cols(); j++)
        {
          // Of course we cannot count the distance of an member to itself
          // since it is =0.0
          if (j == col_queryPt) { continue; }

          // For number of dimensions, add the squared distance of the queryPt
          // to the current point ("j") of each dimensions which equals a
          // squared distance in euclidean space
          temp_distance = 0.0;
          for (size_t l = 0; l < dimension_in; l++)
          {
            temp_distance += pow((input)(row_queryPts[l], j) - (input)(row_queryPts[l], col_queryPt), 2);
          }

          // If we are searching for the actual, "first" nearest neighbor
          // we will compare the distance in "hold_distance", initialized as a huge number,
          // to the calcualted temp_distance ("Is this point nearer to the query point
          // than the previously established nearest point?"). We will, of course,
          // keep the smaller (squared) distance of the two.
          if (i == 0 && temp_distance < hold_distance)
          {
            hold_distance = temp_distance;
          }

          // If we are searching for the "(i + 1)-th" nearest neighbot
          // we will compare "hold_distance" to "temp_distance". If temp_distance is smaller, it
          // will be kept, however, only if it is also larger than the previously established
          // "i-th" nearest neighbor. Otherwise we would always just get the absolute nearest neighbor in
          // a set of points, and not the k-th nearest neighbor.
          else if (i > 0 && temp_distance < hold_distance && temp_distance > distanceList[i - 1])
          {
            hold_distance = temp_distance;
          }
        }
        distanceList[i] = hold_distance;
      }
      float_type keeper = distanceList[k_in - 1u];
      if (buffer == nullptr) delete[] distanceList;
      return keeper;
    }

    float_type knapp_wrapper(Matrix_Class const& input_)
    {
#ifdef _OPENMP
      Matrix_Class input(input_);
#else
      Matrix_Class const& input(input_);
#endif // _OPENMP

      std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Knapp et. al. with corrections (Genome Inform. 2007;18:192-205.)" << std::endl;
      Matrix_Class cov_matr = Matrix_Class{ transposed(input) };
      cov_matr = cov_matr - Matrix_Class( input.cols(), input.cols(), 1. ) * cov_matr / static_cast<float_type>(input.cols());
      cov_matr = transposed(cov_matr) * cov_matr;
      cov_matr *= (1.f / static_cast<float_type>( input.cols() ));
      Matrix_Class eigenvalues;
      Matrix_Class eigenvectors;
	  float_type cov_determ = 0.;
	  int *cov_rank = new int;
      cov_matr.eigensym(eigenvalues, eigenvectors, cov_rank);

      //Remove Eigenvalues that should be zero if cov_matr is singular
      if ((*cov_rank < (int) eigenvalues.rows()) || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90) )
      {
        std::cout << "Notice: covariance matrix is singular, attempting to fix by truncation of Eigenvalues.\n";
        std::cout << "Details: rank of covariance matrix is " << *cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << ".\n";
        if (Config::get().entropy.entropy_remove_dof)
        {
          size_t temp = std::max(6, int((cov_matr.rows() - *cov_rank)));
		      eigenvalues.shed_rows(eigenvalues.rows() - temp, eigenvalues.rows() - 1u);
		      eigenvectors.shed_cols(eigenvectors.cols() - temp, eigenvectors.cols() - 1u);
        }
        else
        {
          eigenvalues.shed_rows((*cov_rank), eigenvalues.rows() - 1u);
          eigenvectors.shed_cols((*cov_rank), eigenvectors.cols() - 1u);
        }
      }
      else if (Config::get().entropy.entropy_remove_dof)
      {
        eigenvectors.shed_cols(0, 5);
        eigenvalues.shed_rows(0, 5);
      }
      delete cov_rank;

      //Calculate PCA Frequencies in quasi-harmonic approximation and Entropy in SHO approximation; provides upper limit of entropy
      Matrix_Class pca_frequencies(eigenvalues.rows(), 1u);
      Matrix_Class alpha_i(pca_frequencies.rows(), 1u);
      Matrix_Class quantum_entropy(pca_frequencies.rows(), 1u);
      float_type entropy_sho = 0;
      for (std::size_t i = 0; i < eigenvalues.rows(); i++)
      {
        pca_frequencies(i, 0u) = sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp / eigenvalues(i, 0u));
        alpha_i(i, 0u) = 1.05457172647 * 10e-34 / (sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp) * sqrt(eigenvalues(i, 0u)));
        quantum_entropy(i, 0u) = ((alpha_i(i, 0u) / (exp(alpha_i(i, 0u)) - 1)) - log(1 - exp(-1 * alpha_i(i, 0u)))) * 1.380648813 * 6.02214129 * 0.239005736;
        entropy_sho += quantum_entropy(i, 0u);
      }
      std::cout << "Entropy in QH-approximation: " << entropy_sho << " cal / (mol * K)\n";
      std::cout << "Starting corrections for anharmonicity and M.I... \n";

      //Corrections for anharmonicity and M.I.
      // I. Create PCA-Modes matrix
      Matrix_Class eigenvectors_t(transposed(eigenvectors));
      Matrix_Class pca_modes = eigenvectors_t * input;
      Matrix_Class entropy_anharmonic(pca_modes.rows(), 1u, 0.);
      Matrix_Class entropy_mi(pca_modes.rows(), pca_modes.rows(), 0.);
      Matrix_Class classical_entropy(pca_modes.rows(), 1u, 0.);
      //std::size_t const size = entropy_anharmonic.rows();


      // II. Calculate Non-Paramteric Entropies
      // Marginal
      const size_t kForKNN = Config::get().entropy.entropy_method_knn_k;
#ifdef _OPENMP
#pragma omp parallel for firstprivate(pca_modes, kForKNN) shared(entropy_anharmonic, entropy_mi)
#endif
      for (int i = 0; i < (int) pca_modes.rows(); i++)
      {
        float_type* buffer = new float_type[kForKNN];
        float_type distance = 0.0;
        {
          for (size_t k = 0; k < pca_modes.cols(); k++)
          {
            distance += log(sqrt(knn_distance(pca_modes, 1, kForKNN, i, k, buffer)));
          }
        }
        distance /= float_type(pca_modes.cols());
        float_type temp = float_type(pca_modes.cols()) * pow(3.14159265358979323846, 0.5);
        temp /= tgamma((1 / 2) + 1);
        distance += log(temp);
        temp = 0;
        if (kForKNN != 1u)
        {
          for (size_t k = 1; k < kForKNN; k++)
          {
            temp += 1.0 / float_type(i);
          }
        }
        distance -= temp;
        distance += 0.5772156649015328606065;
        entropy_anharmonic(i, 0u) = distance;
        distance = 0;

        //MI
        for (size_t j = i + 1; j < pca_modes.rows(); j++)
        {
          std::vector<size_t> query_rows{ (size_t) i,j };
          for (size_t k = 0; k < pca_modes.cols(); k++)
          {
            distance += log(sqrt(knn_distance(pca_modes, 2, kForKNN, query_rows, k, buffer)));
          }

          distance /= 2 * float_type(pca_modes.cols());
          float_type temp2 = float_type(pca_modes.cols()) * pow(3.14159265358979323846, 1);

          distance += log(temp2);
          temp2 = 0;
          if (kForKNN != 1)
          {
            for (size_t u = 1; u < kForKNN; u++)
            {
              temp2 += 1.0 / float_type(u);
            }
          }
          distance -= temp2;
          distance += 0.5772156649015328606065;
          entropy_mi(i, j) = distance;

        }
        delete[] buffer;
      }

      size_t counterForLargeNegativeM_I_Terms = 0u;
      for (size_t i = 0; i < entropy_anharmonic.rows(); i++)
      {
        for (size_t j = (i + 1); j < entropy_anharmonic.rows(); j++)
        {
          if (pca_frequencies(i, 0u) < (Config::get().entropy.entropy_temp * 1.380648813 * 10e-23 / (1.05457172647 * 10e-34)) 
            && pca_frequencies(j, 0u) < (Config::get().entropy.entropy_temp * 1.380648813 * 10e-23 / (1.05457172647 * 10e-34)))
          {
            entropy_mi(i, j) = entropy_anharmonic(i, 0u) + entropy_anharmonic(j, 0u) - entropy_mi(i, j);
            if (entropy_mi(i, j) < 0)
            {
              if (entropy_mi(i, j) < -1)
              {
                counterForLargeNegativeM_I_Terms++;
              }
              entropy_mi(i, j) = 0.0;
            }
            entropy_mi(i, j) *= 1.380648813 * 6.02214129 * 0.239005736;
          }
          else
          {
            if (Config::get().general.verbosity > 4u) std::cout << "Notice: PCA-Modes " << i << " & " << j << " not corrected for M.I. since they are not in the classical limit\n";
            entropy_mi(i, j) = 0.0;
          }
        }
      }
      if (counterForLargeNegativeM_I_Terms > 0u)
      {
        std::cout << "Notice: Large negative M.I. term(s) detected. Check frequency of data sampling. (Do not worry, terms <0.0 are ignored anyway)\n";
      }
      for (size_t i = 0; i < entropy_anharmonic.rows(); i++)
      {
        if (pca_frequencies(i, 0u) < (Config::get().entropy.entropy_temp * 1.380648813 * 10e-23 / (1.05457172647 * 10e-34)))
        {
          classical_entropy(i, 0u) = -1.0 * (log(alpha_i(i, 0u)) + log(sqrt(2 * 3.14159265358979323846 * 2.71828182845904523536)));
          entropy_anharmonic(i, 0u) = classical_entropy(i, 0u) - entropy_anharmonic(i, 0u);
          entropy_anharmonic(i, 0u) *= 1.380648813 * 6.02214129 * 0.239005736;
          if ((entropy_anharmonic(i, 0u) / quantum_entropy(i, 0u)) < 0.007)
          {
            entropy_anharmonic(i, 0u) = 0.0;
            if (Config::get().general.verbosity > 4u) std::cout << "Notice: PCA-Mode " << i << " not corrected for anharmonicity (value too small)";
          }
        }
        else
        {
          if (Config::get().general.verbosity > 4u) std::cout << "Notice: PCA-Mode " << i << " not corrected for anharmonicity since it is not in the classical limit\n";
          entropy_anharmonic(i, 0u) = 0.0;
        }
      }

      // III. Calculate Difference of Entropies
      double delta_entropy = 0;
      for (size_t i = 0; i < entropy_anharmonic.rows(); i++)
      {
        delta_entropy += entropy_anharmonic(i, 0u);
        for (size_t j = (i + 1); j < entropy_anharmonic.rows(); j++)
        {
          delta_entropy += entropy_mi(i, j);
        }
      }
      std::cout << "Correction for entropy: " << delta_entropy << " cal / (mol * K)\n";
      std::cout << "Entropy after correction: " << entropy_sho - delta_entropy << " cal / (mol * K)\n";
      return entropy_sho - delta_entropy;
    }

    float_type hnizdo_wrapper(Matrix_Class const& input_)
    {
#ifdef _OPENMP
      Matrix_Class input(input_);
#else
      Matrix_Class const& input(input_);
#endif // _OPENMP
      std::cout << "\nCommencing entropy calculation:\nNearest-Neighbor Nonparametric Method, according to Hnizdo et al. (DOI: 10.1002/jcc.20589)" << std::endl;
      Matrix_Class marginal_entropy_storage(input.rows(), 1u, 0.);

      const size_t kForKNN = Config::get().entropy.entropy_method_knn_k;
      float_type distance = std::log(std::sqrt(knn_distance(input, input.rows(), kForKNN, (size_t)0u, (size_t)0u)));
#ifdef _OPENMP
#pragma omp parallel firstprivate(input, kForKNN) reduction(+:distance)
      {
#endif
        float_type* buffer = new float_type[kForKNN];
#ifdef _OPENMP
#pragma omp for
#endif
        for (int k = 1; k < (int) input.cols(); k++)
        {
          //Should be equivalent to the original version but with less costly sqrt:
          //This is the squared distance btw
          //distance *= knn_distance(input_workobj,  input.rows(), Config::get().entropy.entropy_method_knn_k, (size_t) 0u, (size_t) k);
          distance += std::log(std::sqrt(knn_distance(input, input.rows(), kForKNN, (size_t)0u, (size_t)k, buffer)));
        }
        delete[] buffer;
#ifdef _OPENMP
      }
#endif

      // Calculate Non-Paramteric Entropies (old)
      /*float_type distance = 0.0;
      for (size_t k = 0; k < this->cols(); k++)
      {
      //Should be equivalent to the original version but with less costly sqrt:
      distance += log(sqrt(this->knn_distance(this->rows(), Config::get().entropy.entropy_method_knn_k, 0, k)));
      }
      distance /= float_type(this->cols());*/

      float_type temp = float_type(input.cols()) * pow(3.14159265358979323846, 0.5);
      temp /= tgamma((1 / 2) + 1);
      distance += log(temp);
      temp = 0;
      if (Config::get().entropy.entropy_method_knn_k != 1)
      {
        for (size_t i = 1; i < Config::get().entropy.entropy_method_knn_k; i++)
        {
          temp += 1.0 / float_type(i);
        }
      }
      distance -= temp;
      distance += 0.5772156649015328606065;

      float_type entropy = distance;
      std::cout << "entropy (dimensionless): " << entropy << '\n';
      std::cout << "entropy: " << entropy * 1.380648813 * 6.02214129 * 0.239005736 << " cal / (mol * K)" << "\n";
      return entropy;
    }

    float_type hnizdo_m_wrapper(Matrix_Class const& input_)
    {
#ifdef _OPENMP
      Matrix_Class input(input_);
#else
      Matrix_Class const& input(input_);
#endif // _OPENMP
      std::cout << "\nCommencing entropy calculation:\nNearest-Neighbor Nonparametric Method - only calculate sum of Marginal Entropies, according to Hnizdo et. al. (DOI: 10.1002/jcc.20589)" << std::endl;
      Matrix_Class marginal_entropy_storage(input.rows(), 1u, 0u);

      //Calculate Non-Paramteric Entropies
      const size_t kForKNN = Config::get().entropy.entropy_method_knn_k;
#ifdef _OPENMP
#pragma omp parallel for firstprivate(input, kForKNN) shared(marginal_entropy_storage)
#endif // _OPENMP
      for (int i = 0; i < (int) input.rows(); i++)
      {
        float_type* buffer = new float_type[kForKNN];
        double distance = 0.0;
        for (size_t k = 0; k < input.cols(); k++)
        {
          distance += log(sqrt(knn_distance(input, 1, kForKNN, (size_t) i, k, buffer))); // set eucledean distance to ouptut
        }
        distance /= float_type(input.cols());

        float_type temp = float_type(input.cols()) * pow(3.14159265358979323846, 0.5);
        temp /= tgamma((1 / 2) + 1);
        distance += log(temp);
        temp = 0;
        if (kForKNN != 1)
        {
          for (size_t k = 1; k < kForKNN; k++)
          {
            temp += 1.0 / float_type(k);
          }
        }
        distance -= temp;
        distance += 0.5772156649015328606065;
        marginal_entropy_storage(i, 0u) = distance;
        delete[] buffer;
      }

      //Calculate sum of of Entropies
      float_type entropy = 0;
      for (size_t i = 0; i < marginal_entropy_storage.rows(); i++)
      {
        entropy += marginal_entropy_storage(i, 0u);
      }
      std::cout << "entropy (dimensionless): " << entropy << '\n';
      std::cout << "entropy: " << entropy * 1.380648813 * 6.02214129 * 0.239005736 << " cal / (mol * K)" << "\n";
      return entropy;
    }

    float_type knapp_m_wrapper(Matrix_Class const& input)
    {
      std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Knapp et. al. without corrections (Genome Inform. 2007;18:192-205.)" << std::endl;
      Matrix_Class cov_matr = (transposed(input));
      cov_matr = cov_matr - Matrix_Class(input.cols(), input.cols(), 1.) * cov_matr / (float_type)input.cols();
      cov_matr = transposed(cov_matr) * cov_matr;
      cov_matr *= (1.f / static_cast<float_type>(input.cols()));
      Matrix_Class eigenvalues;
      Matrix_Class eigenvectors;
	    float_type cov_determ = 0.;
	    int *cov_rank = new int;
      cov_matr.eigensym(eigenvalues, eigenvectors, cov_rank);

      //Remove Eigenvalues that should be zero if cov_matr is singular
      if (( *cov_rank < (int) eigenvalues.rows() ) || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90) )
      {
        std::cout << "Notice: covariance matrix is singular, attempting to fix by truncation of Eigenvalues.\n";
        std::cout << "Details: rank of covariance matrix is " << *cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << ".\n";
        if (Config::get().entropy.entropy_remove_dof)
        {
          size_t temp = std::max(6, int((cov_matr.rows() - *cov_rank)));
          eigenvalues.shed_rows((eigenvalues.rows()) - temp, eigenvalues.rows() - 1u);
          eigenvectors.shed_cols((eigenvectors.cols()) - temp, eigenvectors.cols() - 1u);
        }
        else
        {
		      eigenvalues.shed_rows( (*cov_rank), (eigenvalues.rows()) - 1u);
	    	  eigenvectors.shed_cols( (*cov_rank), (eigenvectors.cols()) - 1u);
        }
      }
      else if (Config::get().entropy.entropy_remove_dof)
      {
        eigenvectors.shed_cols(0, 5);
        eigenvalues.shed_rows(0, 5);
      }
      delete cov_rank;

      //Calculate PCA Frequencies in quasi-harmonic approximation and Entropy in SHO approximation; provides upper limit of entropy
      Matrix_Class pca_frequencies(eigenvalues.rows(), 0u);
      Matrix_Class alpha_i(pca_frequencies.rows(), 0u);
      Matrix_Class quantum_entropy(pca_frequencies.rows(), 0u);
      float_type entropy_sho = 0;
      for (size_t i = 0; i < eigenvalues.rows(); i++)
      {
        pca_frequencies(i, 0u) = sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp / eigenvalues(i, 0u));
        alpha_i(i, 0u) = 1.05457172647 * 10e-34 / (sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp) * sqrt(eigenvalues(i, 0u)));
        quantum_entropy(i, 0u) = ((alpha_i(i, 0u) / (exp(alpha_i(i, 0u)) - 1)) - log(1 - exp(-1 * alpha_i(i, 0u)))) * 1.380648813 * 6.02214129 * 0.239005736;
        entropy_sho += quantum_entropy(i, 0u);
      }
      std::cout << "Entropy in QH-approximation: " << entropy_sho << " cal / (mol * K)\n";
      return entropy_sho;
    }

    float_type karplus_wrapper(Matrix_Class const& input)
    {
      std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Karplus et. al. (DOI 10.1021/ma50003a019)" << std::endl;
      Matrix_Class cov_matr = (transposed(input));
      cov_matr = cov_matr - Matrix_Class( input.cols(), input.cols(), 1. ) * cov_matr / static_cast<float_type>(input.cols());
      cov_matr = transposed(cov_matr) * cov_matr;
      cov_matr = cov_matr / static_cast<float_type>(input.cols());
      float_type entropy = 0.0, cov_determ;
      if (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90)
      {
        std::cout << "Error: Covariance Matrix is singular. Try: a.) using internal coordinates b.) higher sampling rate in MD-Simulation c.) using more advanced methods.\n";
      }
      else if (cov_determ != 0)
      {
        entropy = log(pow(2 * SCON_PI, int(input.rows())) * cov_determ);
        entropy += (input.rows());
        entropy *= 1.380648813 * 6.02214129 * 0.239005736;
        std::cout << "Entropy in classical QH-approximation: " << entropy << " cal / (mol * K)" << std::endl;
      }
      return entropy;
    }

    float_type schlitter_wrapper(Matrix_Class const& input)
    {
      std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Schlitter (see: doi:10.1016/0009-2614(93)89366-P)" << std::endl;
      Matrix_Class cov_matr = transposed(input);
      cov_matr = cov_matr - Matrix_Class(input.cols(), input.cols(), 1.0) * cov_matr / static_cast<float_type>(input.cols());
      cov_matr = transposed(cov_matr) * cov_matr;
      cov_matr = cov_matr / static_cast<float_type>(input.cols());

      cov_matr *= (1.38064813 * /* 10e-23 J/K */ Config::get().entropy.entropy_temp * 2.718281828459 * 2.718281828459 / (1.054571726 /* * 10^-34 Js */ * 1.054571726 * 10e-45));
      cov_matr = cov_matr + Matrix_Class::identity(cov_matr.rows(), cov_matr.cols());
      float_type entropy_sho = cov_matr.determ();

      entropy_sho = log(entropy_sho) * 0.5 * 1.38064813 * 6.02214129 * 0.239;
      //This stems from: k_B * Avogadro * (to_calories) *0.5

      std::cout << "Entropy in QH-approximation: " << entropy_sho << " cal / (mol * K)" << std::endl;
      return entropy_sho;
    }

  }

  /////////////////////////////////////
  //                              /////
  //    E X C L U S I V E L Y     /////
  //      A L I G N M E N T       /////
  //                              /////
  /////////////////////////////////////
  namespace align
  {
    float_type drmsd_calc(coords::Coordinates const& input, coords::Coordinates const& ref)
    {
      if (input.atoms().size() != ref.atoms().size()) throw std::logic_error("Number of atoms of structures passed to drmsd_calc to not match.");

      float_type value = 0;
      for (size_t i = 0; i < input.atoms().size(); i++)
      {
        for (size_t j = 0; j < i; j++)
        {
          float_type holder = sqrt((ref.xyz(i).x() - ref.xyz(j).x()) * (ref.xyz(i).x() - ref.xyz(j).x()) + (ref.xyz(i).y() - ref.xyz(j).y())* (ref.xyz(i).y() - ref.xyz(j).y()) + (ref.xyz(i).z() - ref.xyz(j).z()) * (ref.xyz(i).z() - ref.xyz(j).z()));
          float_type holder2 = sqrt((input.xyz(i).x() - input.xyz(j).x()) * (input.xyz(i).x() - input.xyz(j).x()) + (input.xyz(i).y() - input.xyz(j).y())* (input.xyz(i).y() - input.xyz(j).y()) + (input.xyz(i).z() - input.xyz(j).z()) * (input.xyz(i).z() - input.xyz(j).z()));
          value += (holder2 - holder) * (holder2 - holder);
        }
      }
      return sqrt(value / (double) (input.atoms().size() * (input.atoms().size() + 1u) ) );
    }

    float_type holmsander_calc(coords::Coordinates const& input, coords::Coordinates const& ref, double holmAndSanderDistance)
    {
      if (input.atoms().size() != ref.atoms().size()) throw std::logic_error("Number of atoms of structures passed to drmsd_calc to not match.");

      float_type value(0);
      for (size_t i = 0; i < input.atoms().size(); i++) {
        for (size_t j = 0; j < i; j++)
        {
          float_type holder = sqrt((ref.xyz(i).x() - ref.xyz(j).x()) * (ref.xyz(i).x() - ref.xyz(j).x()) + (ref.xyz(i).y() - ref.xyz(j).y())* (ref.xyz(i).y() - ref.xyz(j).y()) + (ref.xyz(i).z() - ref.xyz(j).z()) * (ref.xyz(i).z() - ref.xyz(j).z()));
          float_type holder2 = sqrt((input.xyz(i).x() - input.xyz(j).x()) * (input.xyz(i).x() - input.xyz(j).x()) + (input.xyz(i).y() - input.xyz(j).y())* (input.xyz(i).y() - input.xyz(j).y()) + (input.xyz(i).z() - input.xyz(j).z()) * (input.xyz(i).z() - input.xyz(j).z()));
          value += abs(holder2 - holder) * exp(-1 * (holder2 + holder)*(holder2 + holder) / (4 * holmAndSanderDistance * holmAndSanderDistance)) / (holder2 + holder);
        }
      }
      return value;
    }

    coords::Coordinates kabschAligned(coords::Coordinates const& inputCoords, coords::Coordinates const& reference, bool centerOfMassAlign)
    {
      coords::Coordinates output(inputCoords);
      if (centerOfMassAlign) centerOfMassAlignment(output);
      kabschAlignment(output, reference);
      return output;
    }

    void kabschAlignment(coords::Coordinates& inputCoords, coords::Coordinates const& reference, bool centerOfMassAlign)
    {
      if (centerOfMassAlign)
      {
        centerOfMassAlignment(inputCoords);
      }

      Matrix_Class input = transfer_to_matr(inputCoords);
      Matrix_Class ref = transfer_to_matr(reference);

      Matrix_Class c(input * transposed(ref));
      //Creates Covariance Matrix

      Matrix_Class s, V, U;
      c.singular_value_decomposition(U, s, V);

      Matrix_Class unit = Matrix_Class::identity(c.rows(), c.rows());
      if ((c.det_sign() < 0)) //Making sure that U will do a proper rotation (rows/columns have to be right handed system)
      {
        unit(2, 2) = -1;
      }
      transpose(U);
      unit = unit * U;
      unit = V * unit;
      input = unit * input;

      inputCoords.set_xyz(transfer_to_3DRepressentation(input));
    }

    void centerOfMassAlignment(coords::Coordinates& coords_in)
    {
      coords::Cartesian_Point com_ref = coords_in.center_of_mass();
      coords_in.move_all_by(-com_ref, true);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                    W R A P P E R F U N C T I O N S                         //
////////////////////////////////////////////////////////////////////////////////

void alignment(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords)
{
  using namespace matop;
  using namespace matop::align;

  coords::Coordinates coordsReferenceStructure(coords), coordsTemporaryStructure(coords);

  // Check if reference structure is in range
  if (Config::get().alignment.reference_frame_num >= ci->size()) throw std::runtime_error("Reference frame number in ALIGN task is bigger than number of frames in the input structure ensemble.");

  auto temporaryPESpoint = ci->PES()[Config::get().alignment.reference_frame_num].structure.cartesian;

  //Alignment to external reference frame (different file)
  if (!Config::get().alignment.align_external_file.empty())
  {
    std::unique_ptr<coords::input::format> externalReferenceStructurePtr(coords::input::new_format());
    coords::Coordinates externalReferenceStructure(externalReferenceStructurePtr->read(Config::get().alignment.align_external_file));
    if (Config::get().alignment.reference_frame_num >= externalReferenceStructurePtr->PES().size())
    {
      throw std::out_of_range("Requested reference frame number not within reference structure ensemble.");
    }
    temporaryPESpoint = externalReferenceStructurePtr->PES()[Config::get().alignment.reference_frame_num].structure.cartesian;
  }
  //Constructs two coordinate objects and sets reference frame according to INPUTFILE
  coordsReferenceStructure.set_xyz(temporaryPESpoint);

  //Construct and Allocate arrays for stringoutput (necessary for OpenMP)
  double mean_value = 0;
  std::string *hold_str, *hold_coords_str;
  hold_str = new std::string[ci->size()];
  hold_coords_str = new std::string[ci->size()];

  //Perform translational alignment for reference frame
  if (Config::get().alignment.traj_align_translational)
  {
    centerOfMassAlignment(coordsReferenceStructure);
  }

  // Output text
  if (Config::get().general.verbosity > 2U) std::cout << "ALIGN preparations done. Starting actual alignment.\n";

#ifdef _OPENMP
  if (Config::get().general.verbosity > 3U) std::cout << "Using openMP for alignment.\n";
  auto const n_omp = static_cast<std::ptrdiff_t>(ci->size());
#pragma omp parallel for firstprivate(coordsReferenceStructure, coordsTemporaryStructure) reduction(+:mean_value) shared(hold_coords_str, hold_str)
  for (std::ptrdiff_t i = 0; i < n_omp; ++i)
#else
  for (std::size_t i = 0; i < ci->size(); ++i)
#endif
  {
    if (i != static_cast<std::ptrdiff_t>(Config::get().alignment.reference_frame_num))
    {
      auto temporaryPESpoint2 = ci->PES()[i].structure.cartesian;
      coordsTemporaryStructure.set_xyz(temporaryPESpoint2);
      //Create temporary objects for current frame

      if (Config::get().alignment.traj_align_translational)
      {
        centerOfMassAlignment(coordsTemporaryStructure);
      }
      if (Config::get().alignment.traj_align_rotational)
      {
        kabschAlignment(coordsTemporaryStructure, coordsReferenceStructure, false);
      }

      if (Config::get().alignment.traj_print_bool)
      {
        if (Config::get().alignment.dist_unit == 0)
          //RMSD
        {
          std::stringstream temporaryStringstream;
          double currentRootMeanSquareDevaition = root_mean_square_deviation(coordsTemporaryStructure.xyz(), coordsReferenceStructure.xyz());
          temporaryStringstream << std::setw(13) << i << " ";
          temporaryStringstream << std::setw(13) << currentRootMeanSquareDevaition << "\n";
          mean_value += currentRootMeanSquareDevaition;
          hold_str[i] = temporaryStringstream.str();
        }
        else if (Config::get().alignment.dist_unit == 1)
          //dRMSD
        {
          std::stringstream temporaryStringstream;
          temporaryStringstream << i << " ";
          double value = (double)drmsd_calc(coordsTemporaryStructure, coordsReferenceStructure);
          temporaryStringstream << std::setw(13) << value << "\n";
          mean_value += value;
          hold_str[i] = temporaryStringstream.str();
        }
        else if (Config::get().alignment.dist_unit == 2)
          //Holm&Sander Distance
        {
          std::stringstream temporaryStringstream;
          double value = (double)holmsander_calc(coordsTemporaryStructure, coordsReferenceStructure, Config::get().alignment.holm_sand_r0);
          temporaryStringstream << std::setw(13) << i << " " << value << "\n";
          mean_value += value;
          hold_str[i] = temporaryStringstream.str();
        }
      }
      //Molecular distance measure calculation

      std::stringstream hold_coords;
      hold_coords << coordsTemporaryStructure;
      hold_coords_str[i] = hold_coords.str();
      //Formatted string-output
    }

    else if (i == static_cast<std::ptrdiff_t>(Config::get().alignment.reference_frame_num))
    {
      std::stringstream hold_coords;
      hold_coords << coordsReferenceStructure;
      hold_coords_str[i] = hold_coords.str();
      //Formatted string-output (first to array because of OpenMP parallelization)
    }
  }

  std::ofstream distance(coords::output::filename("_distances").c_str(), std::ios::app);
  std::ofstream outputstream(coords::output::filename("_aligned").c_str(), std::ios::app);

  if (Config::get().general.verbosity > 2U) std::cout << "Alignment done. Writing structures to file.\n";

  for (size_t i = 0; i < ci->size(); i++)
  {
    if (Config::get().alignment.traj_print_bool)
    {
      distance << hold_str[i];
    }
    outputstream << hold_coords_str[i];
  }
  if (Config::get().alignment.traj_print_bool)
  {
    distance << "\n";
    distance << "Mean value: " << (mean_value / (double)(ci->size() - 1)) << "\n";
  }
  //Formatted string-output

  delete[] hold_str;
  delete[] hold_coords_str;
  //Cleaning Up
}

void entropy(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords)
{
   using namespace matop;
   using namespace matop::entropy;

   //First, adjust number of truncated atoms to be used to zero, in case truncation is not to be used (duh)
   if (!Config::get().entropy.entropy_trunc_atoms_bool) Config::set().entropy.entropy_trunc_atoms_num = std::vector<size_t>();

   //Initialize the reference frame (for alignment etc)
   coords::Coordinates coords_ref(coords);
   auto holder = (*ci).PES()[Config::get().entropy.entropy_ref_frame_num].structure.cartesian;
   coords_ref.set_xyz(holder);

   //Translational alignment of the reference frame
   if (Config::get().entropy.entropy_alignment)
   {
     ::matop::align::centerOfMassAlignment(coords_ref);
   }

   Matrix_Class matrix_aligned;
   //const size_t FRAME_SIZE = int(ci->size());
   if (Config::get().entropy.entropy_alignment && Config::get().entropy.entropy_use_internal)
   {
     std::cout << "Alignment is (in this case) redundant since internal coordinates are used. Alignment is skipped. Check your INPUTFILE please.\n";
     std::cout << "Continuing anyway...";
   }
   //Initializing and checking...

   if (Config::get().entropy.entropy_use_internal)
   {
     matrix_aligned = Matrix_Class((size_t) /* explicitly casting to round down */ ((ci->size() - Config::get().entropy.entropy_start_frame_num) / Config::get().entropy.entropy_offset), \
       Config::get().entropy.entropy_internal_dih.size() * 2u);
   }
   else if (Config::get().entropy.entropy_trunc_atoms_bool)
   {
     matrix_aligned = Matrix_Class((size_t) /* explicitly casting to round down */ ((ci->size() - Config::get().entropy.entropy_start_frame_num) / Config::get().entropy.entropy_offset), \
       Config::get().entropy.entropy_trunc_atoms_num.size() * 3u);
   }
   else
   {
     matrix_aligned = Matrix_Class((size_t) /* explicitly casting to round down */ ((ci->size() - Config::get().entropy.entropy_start_frame_num) / Config::get().entropy.entropy_offset), \
       coords.atoms().size() * 3u);
   }

   /* j counts the (truncated) matrix access, i the frames in ci */
   {
     size_t j = 0;
     if (Config::get().entropy.entropy_use_internal)        //Conversion to internal coordinates if desired
     {
       for (size_t i = Config::get().entropy.entropy_start_frame_num; j < matrix_aligned.rows(); ++j, i += Config::get().entropy.entropy_offset)
       {
         auto holder2 = ci->PES()[i].structure.cartesian;
         coords.set_xyz(holder2);
         matrix_aligned.row(j) = ::matop::transformToOneline(coords, Config::get().entropy.entropy_internal_dih, true);
       }
     }
     else
     {
       for (size_t i = Config::get().entropy.entropy_start_frame_num; j < matrix_aligned.rows(); ++j, i += Config::get().entropy.entropy_offset)
       {
         auto holder2 = ci->PES()[i].structure.cartesian;
         coords.set_xyz(holder2);
         if (Config::get().entropy.entropy_alignment)        //Translational and rotational alignment
         {
           ::matop::align::centerOfMassAlignment(coords); //Alignes center of mass
           ::matop::align::kabschAlignment(coords, coords_ref); //Rotates
         }
         matrix_aligned.row(j) = ::matop::transformToOneline(coords, Config::get().entropy.entropy_trunc_atoms_num, false);
       }
     }
   }

   transpose(matrix_aligned);
   //NECESSARY because of implementation details, don't worry about it for now; rows are DOFs, columns are frames FROM HERE ON!

   if (!Config::get().entropy.entropy_use_internal)
   {
     if (!Config::get().entropy.entropy_trunc_atoms_bool)
     {
       ::matop::massweight(matrix_aligned, coords_ref, true);
     }
     else
     {
       ::matop::massweight(matrix_aligned, coords_ref, true, Config::get().entropy.entropy_trunc_atoms_num);
     }
   }
   //Mass-weightening cartesian coordinates

   for (size_t u = 0u; u < Config::get().entropy.entropy_method.size(); u++)
   {
     Matrix_Class& workobj = matrix_aligned;
     int m = (int)Config::get().entropy.entropy_method[u];
     if (m == 1 || m == 0)
     {
       /*double entropy_value = */karplus_wrapper(workobj);
     }
     if (m == 2)
     {
       /*double entropy_value = */knapp_m_wrapper(workobj);
     }
     if (m == 3 || m == 0)
     {
       /*double entropy_value = */knapp_wrapper(workobj);
     }
     if (m == 4 || m == 0)
     {
       /*double entropy_value = */hnizdo_wrapper(workobj);
     }
     if (m == 5 || m == 0)
     {
       /*double entropy_value = */hnizdo_m_wrapper(workobj);
     }
     if (m == 6 || m == 0)
     {
       /*double entropy_value = */schlitter_wrapper(workobj);
     }
   }
   //Perform desired calculations
}
