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
      coords::Cartesian_Point tempcoord2(input(0, i), input(1, i), input(2, i));
      tempcoord1.push_back(tempcoord2);
    }
    return tempcoord1;
  }

  coords::Representation_Internal transfer_to_internalRepressentation(Matrix_Class const& input)
  {
    coords::Representation_Internal tempcoord1;

    for (size_t i = 0; i < input.cols(); i++)
    {
      coords::internal_type tempcoord2(input(0, i), scon::ang<float_type>(input(1, i)), scon::ang<float_type>(input(2, i)));
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
    if (includedAtoms[includedAtoms.size() - 1] > coords.atoms().size() - 1)
    {
      std::cerr << "You specified a truncation number that is greater than the total number of atoms. Stopping." << std::endl;
      throw;
    }
    if (includedAtoms[0] < 0u)
    {
      std::cerr << "You specified a negative truncation number. Stopping." << std::endl;
      throw;
    }
    else if (internalCoordinates && includedAtoms[0] < 3u)
    {
      std::cerr << "You specified a dihedral with index < 3. Stopping." << std::endl;
      throw;
    }

    if (internalCoordinates)
    {
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
        Matrix_Class transformed_matrix(1u, ((3 * includedAtoms.size())));
        int j = 0;
        size_t quicksearch = 0;
        for (size_t i = 0; i < coords.atoms().size(); i++) //iterate over atoms
        {
          bool checker = false;
          for (size_t l = quicksearch; l < includedAtoms.size(); l++) //iterate over vector of atoms to account for
          {
            if (includedAtoms[l] - 1 == i)
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
            //transformed_matrix(0, j + k) = input(k, i);
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
  //            P C A             /////
  //                              /////
  /////////////////////////////////////
  namespace pca
  {
    void prepare_pca(Matrix_Class const& input, Matrix_Class& eigenvalues, Matrix_Class& eigenvectors, int rank)
    {
      Matrix_Class cov_matr = (transposed(input));
      Matrix_Class ones(input.cols(), input.cols(), 1.0);
      cov_matr = cov_matr - ones * cov_matr / input.cols();
      cov_matr = transposed(cov_matr) * cov_matr;
      cov_matr = cov_matr / input.cols();
      float_type cov_determ = 0.;
      int *cov_rank = new int;
      cov_matr.eigensym(eigenvalues, eigenvectors, cov_rank);
      rank = *cov_rank;
      delete cov_rank;
      if (rank < (int)eigenvalues.rows() || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90))
      {
        std::cout << "Notice: covariance matrix is singular.\n";
        std::cout << "Details: rank of covariance matrix is " << *cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << ".\n";
        if (Config::get().PCA.pca_remove_dof)
        {
          size_t temp = std::max(6, int((cov_matr.rows() - *cov_rank)));
          eigenvalues.shed_rows(eigenvalues.rows() - temp, eigenvalues.rows() - 1);
          eigenvectors.shed_cols(eigenvectors.cols() - temp, eigenvectors.cols() - 1);
        }
      }
      else if (Config::get().PCA.pca_remove_dof)
      {
        eigenvectors.shed_cols(eigenvalues.rows() - 6, eigenvalues.rows() - 1);
        eigenvalues.shed_rows(eigenvectors.cols() - 6, eigenvectors.cols() - 1);
      }
    }

    void output_pca_modes(Matrix_Class& eigenvalues, Matrix_Class& eigenvectors, Matrix_Class& pca_modes, std::string filename, std::string additionalInformation)
    {
      std::cout << "\nWriting PCA-Modes and Eigenvectors of covariance matrix to file\n";

      double sum_of_all_variances = 0.0;
      for (size_t i = 0; i < eigenvalues.rows(); i++)
      {
        sum_of_all_variances += eigenvalues(i);
      }

      if (Config::get().PCA.pca_trunc_var != 1 && Config::get().PCA.pca_trunc_var < 1 && Config::get().PCA.pca_trunc_var > 0)
      {
        double temporary_sum_of_variances = 0.0;
        Matrix_Class pca_modes2 = pca_modes;
        if (eigenvalues.rows() < 2u) throw("Eigenvalues not initialized, error in truncating PCA values. You are in trouble.");
        for (size_t i = 0u - 1; i < eigenvalues.rows(); i++)
        {
          if (Config::get().PCA.pca_trunc_var < temporary_sum_of_variances / sum_of_all_variances)
          {
            size_t rows_of_pca_modes2 = pca_modes2.rows();
            pca_modes2.shed_rows(rows_of_pca_modes2 - i, rows_of_pca_modes2 - 1);
            eigenvalues.shed_rows(rows_of_pca_modes2 - i, rows_of_pca_modes2 - 1);
            break;
          }
          temporary_sum_of_variances += eigenvalues(i);
        }
        pca_modes = pca_modes2;
      }
      else if (Config::get().PCA.pca_trunc_var != 1)
      {
        throw("Error in pca task, check specified trunc_var. (It should be > 0. and < 1.)");
      }
      if (Config::get().PCA.pca_trunc_dim != 0 && Config::get().PCA.pca_trunc_dim < pca_modes.rows())
      {
        Matrix_Class pca_modes2 = pca_modes;
        pca_modes2.shed_rows(Config::get().PCA.pca_trunc_dim, pca_modes.rows() - 1u);
        eigenvalues.shed_rows(Config::get().PCA.pca_trunc_dim, pca_modes.rows() - 1u);
        pca_modes = pca_modes2;
      }
      else if (Config::get().PCA.pca_trunc_dim != 0 && Config::get().PCA.pca_trunc_dim >= pca_modes.rows())
      {
        std::cout << "Notice: Inputfile-specified truncation for PCA is too large, there are less DOFs than the user specified truncation value for dimensionality.\n";
      }
      std::cout << "Working with " << pca_modes.rows() << " dimensions.\n";

      std::ofstream pca_modes_stream(filename, std::ios::out);
      pca_modes_stream << "Working with " << pca_modes.rows() << " dimensions.\n";
      pca_modes_stream << "PCA-Modes are ordered with descending variances.\n\n\n\n";
      pca_modes_stream << "Eigenvalues of Covariance Matrix:\n";
      pca_modes_stream << eigenvalues << "\n\n\n";
      pca_modes_stream << "Eigenvectors of Covariance Matrix:\nSize = " << std::setw(10) <<  eigenvectors.rows() << " x " << std::setw(10) << eigenvectors.cols() << "\n";
      pca_modes_stream << eigenvectors << "\n\n\n";
      pca_modes_stream << "Trajectory in PCA - Modes following (columns are frames, rows are modes)\nSize = " << std::setw(10) << pca_modes.rows() << " x " << std::setw(10) << pca_modes.cols() << "\n";
      pca_modes_stream << pca_modes << "\n";
      pca_modes_stream << "Additional Information following:\n" << additionalInformation << "\n";
    }

    void output_probability_density(Matrix_Class& pca_modes)
    {
      using namespace histo;
      std::vector<size_t> dimensionsToBeUsedInHistogramming = Config::get().PCA.pca_dimensions_for_histogramming;

      if (dimensionsToBeUsedInHistogramming.size() > 2)
      {
        std::cerr << "More than 2 dimensions for histogramming is not yet supported!\n";
        std::cerr << "Working with 2 dimensions";
        dimensionsToBeUsedInHistogramming.resize(2);
      }
      else
      {
        DimensionalHistogram<float_type> *histograms_p;
        if (Config::get().PCA.pca_histogram_number_of_bins > 0u)
        {
          size_t histogramBins = Config::get().PCA.pca_histogram_number_of_bins;
          //histograms_p = new DimensionalHistogram<float_type>((size_t)pca_modes.rows(), histogramBins);
          histograms_p = new DimensionalHistogram<float_type>( (size_t)dimensionsToBeUsedInHistogramming.size(), histogramBins);
        }
        else if (Config::get().PCA.pca_histogram_width > 0.)
        {
          float_type histogramWidth = Config::get().PCA.pca_histogram_width;
          histograms_p = new DimensionalHistogram<float_type>( (size_t)dimensionsToBeUsedInHistogramming.size(), histogramWidth);
        }
        else
        {
          throw "Error in output_probability_density, exiting.\n You need to specify either a numbger of histogram bins or a bin width in the inputfile.\n";
        }

        std::cout << "Starting: Histogramming...";
        //Filling the histogram
        for (size_t j = 0u; j < pca_modes.cols(); j++)
        {
          std::vector<float_type> temp(dimensionsToBeUsedInHistogramming.size());
          for (size_t i = 0u; i < dimensionsToBeUsedInHistogramming.size(); i++)
          {
            temp[i] = pca_modes(dimensionsToBeUsedInHistogramming[i] - 1, j);
          }
          histograms_p->add_value(temp);
        }
        
        histograms_p->distribute();
        histograms_p->writeProbabilityDensity("pca_histograms");
        histograms_p->writeAuxilaryData("pca_histograms_auxdata");
        delete histograms_p;
      }
    }

    void readEigenvectorsAndModes(Matrix_Class& eigenvectors, Matrix_Class& trajectory, std::string& additionalInformation, std::string filename)
    {
      std::ifstream pca_modes_stream(filename, std::ios::in);
      std::string line;
      std::getline(pca_modes_stream, line);
      //int dimensions = std::stoi(line.substr(13, 2));
      while (line.find("Eigenvectors") == std::string::npos)
      {
        std::getline(pca_modes_stream, line);
      }

      std::getline(pca_modes_stream, line);
      eigenvectors.resize(stoi(line.substr(7, 10)), stoi(line.substr(20, 10)));

      for (size_t i = 0u; i < eigenvectors.rows(); i++)
      {
        std::getline(pca_modes_stream, line);
        size_t whitespace = 0u, lastWhitespace = 0u;
        for (size_t j = 0u; j < eigenvectors.cols(); j++)
        {
          lastWhitespace = whitespace;
          whitespace = line.find(" ", lastWhitespace + 1u);

          eigenvectors(i, j) = stod(line.substr(lastWhitespace, whitespace - lastWhitespace));
        }
      }
      std::getline(pca_modes_stream, line);
      std::getline(pca_modes_stream, line);
      std::getline(pca_modes_stream, line);
      std::getline(pca_modes_stream, line);
      std::getline(pca_modes_stream, line);
      trajectory.resize(stoi(line.substr(7, 10)), stoi(line.substr(20, 10)));

      for (size_t i = 0u; i < trajectory.rows(); i++)
      {
        std::getline(pca_modes_stream, line);
        size_t whitespace = 0u, lastWhitespace = 0u;
        for (size_t j = 0u; j < trajectory.cols(); j++)
        {
          lastWhitespace = whitespace;
          whitespace = line.find(" ", lastWhitespace + 1u);

          trajectory(i, j) = stod(line.substr(lastWhitespace, whitespace - lastWhitespace));
        }
      }
      //Additional options following:
      if (std::getline(pca_modes_stream, line))
      {
        std::getline(pca_modes_stream, line);
        if (std::getline(pca_modes_stream, line))
        {
          additionalInformation = line;
        }
        else
        {
          std::cerr << "Could not read additional Information from pca_modes file.\n";
        }
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

    float_type knn_distance(Matrix_Class const& input, size_t const& dimension_in, size_t const& k_in, size_t const& row_queryPt, size_t const& col_queryPt)
    {
      float_type temp_distance = 0.0;
      float_type hold_distance;
      float_type *distanceList = new float_type[k_in];
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
      delete[] distanceList;
      return keeper;
    }

    float_type knn_distance(Matrix_Class const& input, size_t const& dimension_in, size_t const& k_in, std::vector<size_t>& row_queryPts, size_t const& col_queryPt)
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
      float_type *distanceList = new float_type[k_in];
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
      delete[] distanceList;
      return keeper;
    }

    float_type knapp_wrapper(Matrix_Class const& input)
    {
      std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Knapp et. al. with corrections (Genome Inform. 2007;18:192-205.)\n";
      Matrix_Class cov_matr = Matrix_Class{ transposed(input) };
      cov_matr = cov_matr - Matrix_Class( input.cols(), input.cols(), 1. ) * cov_matr / input.cols();
      cov_matr = transposed(cov_matr) * cov_matr;
      cov_matr = cov_matr * (1. / (float_type) input.cols() );
      Matrix_Class eigenvalues;
      Matrix_Class eigenvectors;
	    float_type cov_determ = 0.;
	    int *cov_rank = new int;
      cov_matr.eigensym(eigenvalues, eigenvectors, cov_rank);

      //Remove Eigenvalues that should be zero if cov_matr is singular
      if ((*cov_rank < (int) eigenvalues.rows()) || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90) )
      {
        std::cout << "Notice: covariance matrix is singular, attempting to fix by truncation of Eigenvalues." << lineend;
        std::cout << "Details: rank of covariance matrix is " << *cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << "." << lineend;
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


      //Calculate PCA Frequencies in quasi-harmonic approximation and Entropy in SHO approximation; provides upper limit of entropy
      Matrix_Class pca_frequencies(eigenvalues.rows());
      Matrix_Class alpha_i(pca_frequencies.rows());
      Matrix_Class quantum_entropy(pca_frequencies.rows());
      float_type entropy_sho = 0;
      for (std::size_t i = 0; i < eigenvalues.rows(); i++)
      {
        pca_frequencies(i) = sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp / eigenvalues(i));
        alpha_i(i) = 1.05457172647 * 10e-34 / (sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp) * sqrt(eigenvalues(i)));
        quantum_entropy(i) = ((alpha_i(i) / (exp(alpha_i(i)) - 1)) - log(1 - exp(-1 * alpha_i(i)))) * 1.380648813 * 6.02214129 * 0.239005736;
        entropy_sho += quantum_entropy(i);
      }
      std::cout << "Entropy in QH-approximation: " << entropy_sho << " cal / (mol * K)" << lineend;
      std::cout << "Starting corrections for anharmonicity and M.I... " << lineend;

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
      for (size_t i = 0; i < pca_modes.rows(); i++)
      {
        float_type distance = 0.0;
        for (size_t k = 0; k < pca_modes.cols(); k++)
        {
          distance += log(sqrt(knn_distance(pca_modes, 1, Config::get().entropy.entropy_method_knn_k, i, k)));
        }
        distance /= float_type(pca_modes.cols());
        float_type temp = float_type(pca_modes.cols()) * pow(3.14159265358979323846, 0.5);
        temp /= tgamma((1 / 2) + 1);
        distance += log(temp);
        temp = 0;
        if (Config::get().entropy.entropy_method_knn_k != 1u)
        {
          for (size_t k = 1; k < Config::get().entropy.entropy_method_knn_k; k++)
          {
            temp += 1.0 / float_type(i);
          }
        }
        distance -= temp;
        distance += 0.5772156649015328606065;
        entropy_anharmonic(i) = distance;
        distance = 0;

        //MI
        for (size_t j = i + 1; j < pca_modes.rows(); j++)
        {
          std::vector<size_t> query_rows{ i,j };
          for (size_t k = 0; k < pca_modes.cols(); k++)
          {
            distance += log(sqrt(knn_distance(pca_modes, 2, Config::get().entropy.entropy_method_knn_k, query_rows, k)));
          }

          distance /= 2 * float_type(pca_modes.cols());
          float_type temp2 = float_type(pca_modes.cols()) * pow(3.14159265358979323846, 1);

          distance += log(temp2);
          temp2 = 0;
          if (Config::get().entropy.entropy_method_knn_k != 1)
          {
            for (size_t u = 1; u < Config::get().entropy.entropy_method_knn_k; u++)
            {
              temp2 += 1.0 / float_type(u);
            }
          }
          distance -= temp2;
          distance += 0.5772156649015328606065;
          entropy_mi(i, j) = distance;

        }
      }

      size_t counterForLargeNegativeM_I_Terms = 0u;
      for (size_t i = 0; i < entropy_anharmonic.rows(); i++)
      {
        for (size_t j = (i + 1); j < entropy_anharmonic.rows(); j++)
        {
          if (pca_frequencies(i) < (Config::get().entropy.entropy_temp * 1.380648813 * 10e-23 / (1.05457172647 * 10e-34)) && pca_frequencies(j) < (Config::get().entropy.entropy_temp * 1.380648813 * 10e-23 / (1.05457172647 * 10e-34)))
          {
            entropy_mi(i, j) = entropy_anharmonic(i) + entropy_anharmonic(j) - entropy_mi(i, j);
            if (entropy_mi(i, j) < 0)
            {
              if (entropy_mi(i, j) < -1)
              {
                counterForLargeNegativeM_I_Terms++;
                //std::cout << "Notice: Large negative M.I. term detected (value: " << entropy_mi(i, j) << "). Check frequency of data sampling." << lineend;
              }
              entropy_mi(i, j) = 0.0;
            }
            entropy_mi(i, j) *= 1.380648813 * 6.02214129 * 0.239005736;
          }
          else
          {
            std::cout << "Notice: PCA-Modes " << i << " & " << j << " not corrected for M.I. since they are not in the classical limit" << lineend;
            entropy_mi(i, j) = 0.0;
          }
        }
      }
      if (counterForLargeNegativeM_I_Terms > 0u)
      {
        std::cout << "Notice: Large negative M.I. term(s) detected. Check frequency of data sampling. (Do not worry, terms <0.0 are ignored anyway)" << lineend;
      }
      for (size_t i = 0; i < entropy_anharmonic.rows(); i++)
      {
        if (pca_frequencies(i) < (Config::get().entropy.entropy_temp * 1.380648813 * 10e-23 / (1.05457172647 * 10e-34)))
        {
          classical_entropy(i) = -1.0 * (log(alpha_i(i)) + log(sqrt(2 * 3.14159265358979323846 * 2.71828182845904523536)));
          entropy_anharmonic(i) = classical_entropy(i) - entropy_anharmonic(i);
          entropy_anharmonic(i) *= 1.380648813 * 6.02214129 * 0.239005736;
          if ((entropy_anharmonic(i) / quantum_entropy(i)) < 0.007)
          {
            entropy_anharmonic(i) = 0.0;
            std::cout << "Notice: PCA-Mode " << i << " not corrected for anharmonicity (value too small)";
          }
        }
        else
        {
          std::cout << "Notice: PCA-Mode " << i << " not corrected for anharmonicity since it is not in the classical limit" << lineend;
          entropy_anharmonic(i) = 0.0;
        }
      }

      // III. Calculate Difference of Entropies
      double delta_entropy = 0;
      for (size_t i = 0; i < entropy_anharmonic.rows(); i++)
      {
        delta_entropy += entropy_anharmonic(i);
        for (size_t j = (i + 1); j < entropy_anharmonic.rows(); j++)
        {
          delta_entropy += entropy_mi(i, j);
        }
      }
      std::cout << "Correction for entropy: " << delta_entropy << " cal / (mol * K)" << lineend;
      std::cout << "Entropy after correction: " << entropy_sho - delta_entropy << " cal / (mol * K)" << lineend;
      return entropy_sho - delta_entropy;
    }

    float_type hnizdo_wrapper(Matrix_Class const& input)
    {
      std::cout << "\nCommencing entropy calculation:\nNearest-Neighbor Nonparametric Method, according to Hnizdo et al. (DOI: 10.1002/jcc.20589)\n";
      Matrix_Class marginal_entropy_storage(input.rows(), 1u, 0.);

      float_type distance = knn_distance(input, input.rows(), Config::get().entropy.entropy_method_knn_k, (size_t) 0u, (size_t) 1u);
      for (size_t k = 1; k < input.cols(); k++)
      {
        //Should be equivalent to the original version but with less costly sqrt:
        distance *= knn_distance(input,  input.rows(), Config::get().entropy.entropy_method_knn_k, (size_t) 0u, (size_t) k);
      }
      distance = std::log(std::sqrt(distance));

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
      std::cout << "entropy (dimensionless): " << entropy << lineend;
      std::cout << "entropy: " << entropy * 1.380648813 * 6.02214129 * 0.239005736 << " cal / (mol * K)" << "\n";
      return entropy;
    }

    float_type hnizdo_m_wrapper(Matrix_Class const& input)
    {
      std::cout << "\nCommencing entropy calculation:\nNearest-Neighbor Nonparametric Method - only calculate sum of Marginal Entropies, according to Hnizdo et. al. (DOI: 10.1002/jcc.20589)\n";
      Matrix_Class marginal_entropy_storage(input.rows(), 1u, 0u);

      //Calculate Non-Paramteric Entropies
      for (size_t i = 0; i < input.rows(); i++)
      {
        double distance = 0.0;
        for (size_t k = 0; k < input.cols(); k++)
        {
          distance += log(sqrt(knn_distance(input, 1, Config::get().entropy.entropy_method_knn_k, i, k))); // set eucledean distance to ouptut
        }
        distance /= float_type(input.cols());

        float_type temp = float_type(input.cols()) * pow(3.14159265358979323846, 0.5);
        temp /= tgamma((1 / 2) + 1);
        distance += log(temp);
        temp = 0;
        if (Config::get().entropy.entropy_method_knn_k != 1)
        {
          for (size_t k = 1; k < Config::get().entropy.entropy_method_knn_k; k++)
          {
            temp += 1.0 / float_type(k);
          }
        }
        distance -= temp;
        distance += 0.5772156649015328606065;
        marginal_entropy_storage(i) = distance;
      }

      //Calculate Difference of Entropies
      float_type entropy = 0;
      for (size_t i = 0; i < marginal_entropy_storage.rows(); i++)
      {
        entropy += marginal_entropy_storage(i);
      }
      std::cout << "entropy (dimensionless): " << entropy << lineend;
      std::cout << "entropy: " << entropy * 1.380648813 * 6.02214129 * 0.239005736 << " cal / (mol * K)" << "\n";
      return entropy;
    }

    float_type knapp_m_wrapper(Matrix_Class const& input)
    {
      std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Knapp et. al. without corrections (Genome Inform. 2007;18:192-205.)\n";
      Matrix_Class cov_matr = (transposed(input));
      cov_matr = cov_matr - Matrix_Class(input.cols(), input.cols(), 1.) * cov_matr / (float_type)input.cols();
      cov_matr = transposed(cov_matr) * cov_matr;
      cov_matr = cov_matr * (1.0 / (float_type)input.cols());
      Matrix_Class eigenvalues;
      Matrix_Class eigenvectors;
	    float_type cov_determ = 0.;
	    int *cov_rank = new int;
      cov_matr.eigensym(eigenvalues, eigenvectors, cov_rank);

      //Remove Eigenvalues that should be zero if cov_matr is singular
      if (( *cov_rank < (int) eigenvalues.rows() ) || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90) )
      {
        std::cout << "Notice: covariance matrix is singular, attempting to fix by truncation of Eigenvalues." << lineend;
        std::cout << "Details: rank of covariance matrix is " << *cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << "." << lineend;
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


      //Calculate PCA Frequencies in quasi-harmonic approximation and Entropy in SHO approximation; provides upper limit of entropy
      Matrix_Class pca_frequencies(eigenvalues.rows());
      Matrix_Class alpha_i(pca_frequencies.rows());
      Matrix_Class quantum_entropy(pca_frequencies.rows());
      float_type entropy_sho = 0;
      for (size_t i = 0; i < eigenvalues.rows(); i++)
      {
        pca_frequencies(i) = sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp / eigenvalues(i));
        alpha_i(i) = 1.05457172647 * 10e-34 / (sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp) * sqrt(eigenvalues(i)));
        quantum_entropy(i) = ((alpha_i(i) / (exp(alpha_i(i)) - 1)) - log(1 - exp(-1 * alpha_i(i)))) * 1.380648813 * 6.02214129 * 0.239005736;
        entropy_sho += quantum_entropy(i);
      }
      std::cout << "Entropy in QH-approximation: " << entropy_sho << " cal / (mol * K)" << lineend;
      return entropy_sho;
    }

    float_type karplus_wrapper(Matrix_Class const& input)
    {
      std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Karplus et. al. (DOI 10.1021/ma50003a019)\n";
      Matrix_Class cov_matr = (transposed(input));
      cov_matr = cov_matr - Matrix_Class( input.cols(), input.rows(), 1. ) * cov_matr / input.cols();
      cov_matr = transposed(cov_matr) * cov_matr;
      cov_matr = cov_matr / input.cols();
      float_type entropy = 0.0, cov_determ;
      if (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90)
      {
        std::cout << "Error: Covariance Matrix is singular. Try: a.) using internal coordinates b.) higher sampling rate in MD-Simulation c.) using more advanced methods.\n";
      }
      else if (cov_determ != 0)
      {
        entropy = log(pow(2 * PI, int(input.rows())) * cov_determ);
        entropy += (input.rows());
        entropy *= 1.380648813 * 6.02214129 * 0.239005736;
        std::cout << "Entropy in classical QH-approximation: " << entropy << " cal / (mol * K)" << std::endl;
      }
      return entropy;
    }

    float_type schlitter_wrapper(Matrix_Class const& input)
    {
      std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Schlitter (see: doi:10.1016/0009-2614(93)89366-P)\n";
      Matrix_Class cov_matr = transposed(input);
      cov_matr = cov_matr - Matrix_Class(input.cols(), input.cols(), 1.0) * cov_matr / input.cols();
      cov_matr = transposed(cov_matr) * cov_matr;
      cov_matr = cov_matr / input.cols();

      /*
      //Alternative method for calculation, is mathematically equivalent
      arma::Col<float_type> arma_eigenvalues = eig_sym(arma_cov_matr);
      float_type entropy_new = 1.0;
      for (int i = 0; i < arma_eigenvalues.n_rows; i++)
      {
      std::cout << arma_eigenvalues(i) << "\n";
      arma_eigenvalues(i) *= 1.38064813 * Config::get().entropy.entropy_temp * 2.718281828459 * 2.718281828459 / (1.054571726 * 1.054571726 * 10e-45);
      std::cout << arma_eigenvalues(i) << "\n";
      entropy_new *= 1 + arma_eigenvalues(i);
      }
      std::cout << entropy_new << "\n";
      entropy_new = log(entropy_new);
      std::cout << entropy_new << "\n";
      entropy_new = 0.5 * 1.38064813 * 10e-23;
      std::cout << entropy_new << "\n";
      */

      cov_matr = cov_matr * (1.38064813 * Config::get().entropy.entropy_temp * 2.718281828459 * 2.718281828459 / (1.0546 * 1.0546 * 10e-45));
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
  if (Config::get().alignment.traj_align_bool)
  {
    centerOfMassAlignment(coordsReferenceStructure);
  }

#ifdef _OPENMP
#pragma omp parallel for firstprivate(coords, coordsReferenceStructure, coordsTemporaryStructure, matrixReferenceStructure) reduction(+:mean_value) shared(hold_coords_str, hold_str)
#endif
  for (int i = 0; i < (int)ci->size(); i++)
  {
    if (i != Config::get().alignment.reference_frame_num)
    {
      auto temporaryPESpoint2 = ci->PES()[i].structure.cartesian;
      coordsTemporaryStructure.set_xyz(temporaryPESpoint2);
      //Create temporary objects for current frame

      if (Config::get().alignment.traj_align_bool)
      {
        kabschAlignment(coordsTemporaryStructure, coordsReferenceStructure, true);
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

    else if (i == Config::get().alignment.reference_frame_num)
    {
      std::stringstream hold_coords;
      hold_coords << coordsReferenceStructure;
      hold_coords_str[i] = hold_coords.str();
      //Formatted string-output (first to array because of OpenMP parallelization)
    }
  }

  std::ofstream distance(coords::output::filename("_distances").c_str(), std::ios::app);
  std::ofstream outputstream(coords::output::filename("_aligned").c_str(), std::ios::app);
  for (size_t i = 0; i < ci->size(); i++)
  {
    if (Config::get().alignment.traj_print_bool)
    {
      distance << hold_str[i];
    }
    if (Config::get().alignment.traj_align_bool)
    {
      outputstream << hold_coords_str[i];
    }
  }
  distance << "\n";
  distance << "Mean value: " << (mean_value / (double) (ci->size() - 1)) << "\n";
  //Formatted string-output

  delete[] hold_str;
  delete[] hold_coords_str;
  //Cleaning Up
}

void pca_gen(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords)
{
  using namespace matop;
  using namespace matop::pca;

  Matrix_Class pca_modes, eigenvalues, eigenvectors, matrix_aligned;

  if (!Config::get().PCA.pca_read_modes)
  {
    coords::Coordinates coords_ref(coords);
    auto holder = ci->PES()[Config::get().PCA.pca_ref_frame_num].structure.cartesian;
    coords_ref.set_xyz(holder);
    //Constructs two coordinate objects and sets reference frame according to INPUTFILE

    //Perform translational alignment for reference frame
    if (Config::get().PCA.pca_alignment)
    {
      align::centerOfMassAlignment(coords_ref);
    }

    //bool has_it_started = false;
    const size_t FRAME_SIZE = (size_t)ci->size();
    //Initializing some stuff

    for (size_t i = 0; i < FRAME_SIZE; i++)
    {
      if ((Config::get().PCA.pca_start_frame_num <= i) && ((Config::get().PCA.pca_start_frame_num % Config::get().PCA.pca_offset) == (i % Config::get().PCA.pca_offset)))
      {
        auto holder2 = ci->PES()[i].structure.cartesian;
        coords.set_xyz(holder2);
        //Initializing current frame

        if (Config::get().PCA.pca_alignment && !Config::get().PCA.pca_use_internal)
        {
          align::centerOfMassAlignment(coords); //Alignes center of mass
          align::kabschAlignment(coords, coords_ref); //Rotates
        }
        //Translational and rotational alignment

        if (Config::get().PCA.pca_use_internal)
        {
          coords.to_internal();
        }
        //Conversion to internal coordinates if desired

        if (Config::get().PCA.pca_start_frame_num < i)
        {
          if (Config::get().PCA.pca_use_internal)
          {
            matrix_aligned.append_bottom(transformToOneline(coords, Config::get().PCA.pca_internal_dih, true));
          }
          else
          {
            //Works for full and truncated PCA
            matrix_aligned.append_bottom(transformToOneline(coords, Config::get().PCA.pca_trunc_atoms_num, false));
          }
        }
        else if (Config::get().PCA.pca_start_frame_num == i)
        {
          if (Config::get().PCA.pca_use_internal)
          {
            //First a little check if the user-specified values are reasonable
            matrix_aligned = transformToOneline(coords, Config::get().PCA.pca_internal_dih, true);
          }
          else
          {
            matrix_aligned = transformToOneline(coords, Config::get().PCA.pca_trunc_atoms_num, false);
          }
        }
        //Building one huge [frames] x [coordinates] matrix by appending for every frame
      }
    }

    transpose(matrix_aligned);
    // NECESSARY because of implementation details, don't worry about it for now
    // Matrix is now [coordinates] x [frames] !!! !!!
    // THIS IS THE CONVENTION FOR USAGE, STICK TO IT!

    if (!Config::get().PCA.pca_use_internal)
    {
      coords_ref.set_xyz(holder);
      if (!Config::get().PCA.pca_trunc_atoms_bool)
      {
        massweight(matrix_aligned, coords_ref, false);
      }
      else
      {
        massweight(matrix_aligned, coords_ref, false, Config::get().PCA.pca_trunc_atoms_num);
      }
    }
    //Mass-weightening coordinates if cartesians are used
    prepare_pca(matrix_aligned, eigenvalues, eigenvectors);
  }

  if (Config::get().PCA.pca_read_vectors)
  {
    std::string iAmNotImportant_YouMayDiscardMe;
    readEigenvectorsAndModes(eigenvectors, pca_modes, iAmNotImportant_YouMayDiscardMe);
  }

  if (!Config::get().PCA.pca_read_modes)
  {
    pca_modes = transposed(eigenvectors) * matrix_aligned;
  }

  ///////////////////////////////////////

  std::string additionalInformation;
  if (Config::get().PCA.pca_use_internal)
  {
    additionalInformation += "int ";
    for (size_t i = 0u; i < Config::get().PCA.pca_internal_dih.size(); i++)
    {
      additionalInformation += "d" + std::to_string(Config::get().PCA.pca_internal_dih[i]) + " ";
    }
  }
  else
  {
    additionalInformation += "car ";
    for (size_t i = 0u; Config::get().PCA.pca_trunc_atoms_bool && i < Config::get().PCA.pca_trunc_atoms_num.size(); i++)
    {
      additionalInformation += std::to_string(Config::get().PCA.pca_trunc_atoms_num[i]) + " ";
    }
  }

  ///////////////////////////////////////

  if (Config::get().PCA.pca_read_vectors)
  {
    output_pca_modes(eigenvalues, eigenvectors, pca_modes, "pca_modes_new.dat", additionalInformation);
  }
  else
  {
    output_pca_modes(eigenvalues, eigenvectors, pca_modes, "pca_modes.dat", additionalInformation);
  }

  if (Config::get().PCA.pca_print_probability_density)
  {
    output_probability_density(pca_modes);
  }

  ///////////////////////////////////////
}

void pca_proc(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords)
{
  using namespace matop;
  using namespace matop::pca;
  Matrix_Class eigenvectors, trajectory;
  std::vector<size_t> structuresToBeWrittenToFile;
  std::string additionalInformation;
  readEigenvectorsAndModes(eigenvectors, trajectory, additionalInformation);
  if (Config::get().PCA.proc_desired_start.size() > trajectory.rows() || Config::get().PCA.proc_desired_stop.size() > trajectory.rows())
  {
    std::cerr << "Desired PCA-Ranges have higher dimensionality then modes. Omitting the last values.\n";
  }

  for (size_t j = 0u; j < trajectory.cols(); j++)
  {
    bool isStructureInRange = true;
    for (size_t i = 0u; i < trajectory.rows() && i < std::max(Config::get().PCA.proc_desired_stop.size(), Config::get().PCA.proc_desired_start.size()); i++)

    {
      if (i < Config::get().PCA.proc_desired_start.size())
      {
        if (trajectory(i, j) < Config::get().PCA.proc_desired_start[i])
        {
          isStructureInRange = false;
        }
      }
      if (i < Config::get().PCA.proc_desired_stop.size())
      {
        if (trajectory(i, j) > Config::get().PCA.proc_desired_stop[i])
        {
          isStructureInRange = false;
        }
      }
    }
    if (isStructureInRange)
    {
      structuresToBeWrittenToFile.push_back(j);
    }
  }

  if (Config::get().general.verbosity >= 3u) std::cout << "Found " << structuresToBeWrittenToFile.size() << " structures in desired range.\n";

  //Undoing PCA
  trajectory = eigenvectors * trajectory;

  std::ofstream outstream(coords::output::filename("_pca_selection").c_str(), std::ios::app);

  // Case: Cartesian Coordinates.
  if (additionalInformation.substr(0, 3) == "car")
  {
    //Additional Information Processing -> Read "DOFS that were used" from file. Put their identifying numbers in vector.
    std::stringstream ss(additionalInformation.substr(4));
    std::string buffer;
    std::vector<size_t> tokens;
    std::deque<bool> alreadyFoundStructures(ci->size(), false);

    while (ss >> buffer) tokens.push_back((size_t)std::stoi(buffer));

    if (tokens.size() != 0u)
    {
      undoMassweight(trajectory, coords, false, tokens);
      for (size_t i = 0u; i < structuresToBeWrittenToFile.size(); i++)
      {
        Matrix_Class out_mat(3, trajectory.rows() / 3u);
        for (size_t j = 0u; j < trajectory.rows(); j = j + 3)
        {
          out_mat(0, j / 3u) = trajectory(j, structuresToBeWrittenToFile[i]);
          out_mat(1, j / 3u) = trajectory(j + 1u, structuresToBeWrittenToFile[i]);
          out_mat(2, j / 3u) = trajectory(j + 2u, structuresToBeWrittenToFile[i]);
        }
        coords::Coordinates current(coords);

        // For every partial (truncated) structure that is inside the user-defined 
        // range regarding its PCA-Modes,
        // we now search the matching full structure in the input trajectory.
        bool structureFound = false;
        auto structureCartesian = ci->PES()[0].structure.cartesian;
        int structureNumber = -1;

        for (size_t k = 0u; k < ci->size() && !structureFound; k++)
        {
          if (alreadyFoundStructures[k]) continue;

          structureCartesian = ci->PES()[k].structure.cartesian; //Current structure
          structureFound = true;
          structureNumber = (int)k;
          // Remeber, we are now iterating over certain atoms (those to which
          // the PCA was truncated.
          for (size_t l = 0u; l < tokens.size(); l++)
          {
            // If abs() of diff of every coordinate is smaller than 0.5% of coordinate (or, if this value
            // is very small, the arbitrary cutoff 2e-4), consider it a match.
            // However, we look for "not-matching" and break the loop. If everything matches, we continue. 
            // Thats why we negate the criterion in the if clause (!)
            float_type xCompare = 0.005 * std::max(std::abs(structureCartesian[tokens[l]].x()), 2e-4);
            float_type yCompare = 0.005 * std::max(std::abs(structureCartesian[tokens[l]].y()), 2e-4);
            float_type zCompare = 0.005 * std::max(std::abs(structureCartesian[tokens[l]].z()), 2e-4);

            if (!(std::abs(out_mat(0, l) - structureCartesian[tokens[l]].x()) <= xCompare &&
              std::abs(out_mat(1, l) - structureCartesian[tokens[l]].y()) <= yCompare &&
              std::abs(out_mat(2, l) - structureCartesian[tokens[l]].z()) <= zCompare))
            {
              structureFound = false;
              break;
            }
          }
        }
        if (structureFound)
        {
          // Horray, we found it, now write it out!
          alreadyFoundStructures[structureNumber] = true;
          current.set_xyz(structureCartesian);
          current.to_internal();
          outstream << current;
        }
        else
        {
          std::cerr << "Could not find structure restored from PCA-Modes in ensemble of structures from original coordinates.\n";
          std::cerr << "This means that there was no provided structure with a deviance of less than 0.5% to the current restored structure.\n\n";
        }
      }
    }
    else
    {
      // Here, we merely restore the coordinates from the PCA-modes
      // since no trucnation took plage, and write them out.
      undoMassweight(trajectory, coords, false);
      for (size_t i = 0u; i < structuresToBeWrittenToFile.size(); i++)
      {
        Matrix_Class out_mat(3, trajectory.rows() / 3u);
        for (size_t j = 0u; j < trajectory.rows(); j = j + 3)
        {
          out_mat(0, j / 3u) = trajectory(j, structuresToBeWrittenToFile[i]);
          out_mat(1, j / 3u) = trajectory(j + 1u, structuresToBeWrittenToFile[i]);
          out_mat(2, j / 3u) = trajectory(j + 2u, structuresToBeWrittenToFile[i]);
        }
        coords::Coordinates out(coords);
        out.set_xyz(transfer_to_3DRepressentation(out_mat));
        out.to_internal();
        outstream << out;
      }
    }
  }

  // Case: Internal Coordiantes
  else if (additionalInformation.substr(0, 3) == "int")
  {
    // [0]: bond distance tokens, [1]: bond angle tokens [2]: dihedrals tokens.
    // Here we keep track of which DOFs are to be considered
    std::deque<bool> tokens(coords.atoms().size(), false);
    std::deque<bool> alreadyFoundStructures(ci->size(), false);

    //Additional Information Processing -> Read "DOFS that were used" from file. Store in "tokens".
    std::stringstream ss(additionalInformation.substr(4));
    std::string buffer;
    // Read the additional information
    while (ss >> buffer)
    {
      if (buffer.substr(0, 1) == "d") tokens[(size_t)std::stoi(buffer.substr(1))] = true;
    }

    // This is gonna be complicated. Im sorry.
    // Iterate over structures chosen from PCA-Ensemble
    for (size_t i = 0u; i < structuresToBeWrittenToFile.size(); i++)
    {
      bool structureFound = false;
      //auto structureCartesian = ci->PES()[0].structure.cartesian;
      int structureNumber = -1;

      //Iterate over input strucutures -> find matching structure
      for (size_t k = 0u; k < ci->size() && !structureFound; k++)
      {
        if (alreadyFoundStructures[k]) continue;
        coords.set_internal(ci->PES()[k].structure.intern);
        size_t quicksearch = 0u;
        structureFound = true;
        structureNumber = (int)k;
        //Iterating over atoms, see if they all match
        for (size_t j = 0u; j < tokens.size(); j++)
        {
          // If abs() of diff of every coordinate is smaller than 0.1% of coordinate, consider it a match.
          // However, we look for "not-matching" and break the loop. If everything matches, we continue. 
          // Thats why we negate the criterion in the if clause (!)
          if (tokens[j] == true)
          {
            auto compareFromPCA1 = trajectory(quicksearch, structuresToBeWrittenToFile[i]);
            auto compareFromTrajectory1 = std::cos(coords.intern(j).azimuth().radians());
            auto compareFromPCA2 = trajectory(quicksearch + 1u, structuresToBeWrittenToFile[i]);
            auto compareFromTrajectory2 = std::sin(coords.intern(j).azimuth().radians());
            bool found1 = std::abs(compareFromTrajectory1 - compareFromPCA1) <= 0.001 * std::abs(compareFromPCA1) \
              || std::abs(compareFromTrajectory1 - compareFromPCA1) < 0.0000001;
            bool found2 = std::abs(compareFromTrajectory2 - compareFromPCA2) <= 0.001 * std::abs(compareFromPCA2) \
              || std::abs(compareFromTrajectory2 - compareFromPCA2) < 0.0000001;
            //If structures did not match
            if (!(found1 && found2))
            {
              structureFound = false;
              break;
            }
            quicksearch += 2u;
          }
        }
        //If match was found, write out.
        if (structureFound)
        {
          alreadyFoundStructures[structureNumber] = true;
          coords.set_pes(ci->PES()[k]);
          outstream << coords;
        }
      }
      if (!structureFound)
      {
        std::cerr << "Could not find structure restored from PCA-Modes in ensemble of structures from original coordinates.\n";
        std::cerr << "You probably made a mistake somewhere in your INPUTFILE.\nDo not consider structures written out after this message as valid.\n";
      }
    }
  }
}

void entropy(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords)
{
   using namespace matop;
   using namespace matop::entropy;

   //Initialize the reference frame (for alignment etc)
   coords::Coordinates coords_ref(coords);
   auto holder = (*ci).PES()[Config::get().entropy.entropy_ref_frame_num].structure.cartesian;
   coords_ref.set_xyz(holder);

   //Translational alignment of the reference frame
   if (Config::get().entropy.entropy_alignment)
   {
     align::centerOfMassAlignment(coords_ref);
   }

   Matrix_Class matrix_aligned;
   const size_t FRAME_SIZE = int(ci->size());
   if (Config::get().entropy.entropy_alignment && Config::get().entropy.entropy_use_internal)
   {
     std::cerr << "Alignment is (in this case) redundant since internal coordinates are used. Alignment is skipped. Check your INPUTFILE please.\n";
     std::cerr << "Continuing anyway...";
   }
   //Initializing and checking...

   for (size_t i = 0; i < FRAME_SIZE; i++)
   {
     // Meet-your-Maker-Note: Keep it like this with the seeminlgy stupid "is it started", please.
     if (Config::get().entropy.entropy_start_frame_num <= i && (Config::get().entropy.entropy_start_frame_num % Config::get().entropy.entropy_offset == i % Config::get().entropy.entropy_offset))
     {
       auto holder2 = ci->PES()[i].structure.cartesian;
       coords.set_xyz(holder2);
       //Initializing current frame

       if (Config::get().entropy.entropy_alignment && !Config::get().entropy.entropy_use_internal)
       {
         align::centerOfMassAlignment(coords); //Alignes center of mass
         align::kabschAlignment(coords, coords_ref); //Rotates
       }
       //Translational and rotational alignment

       if (Config::get().entropy.entropy_use_internal)
       {
         coords.to_internal();
       }
       //Conversion to internal coordinates if desired

       if (Config::get().entropy.entropy_start_frame_num < i)
       {
         if (Config::get().entropy.entropy_use_internal)
         {
           matrix_aligned.append_bottom(transformToOneline(coords, Config::get().entropy.entropy_internal_dih, true));
         }
         else
         {
           matrix_aligned.append_bottom(transformToOneline(coords, Config::get().entropy.entropy_trunc_atoms_num, false));
         }
       }
       else if (Config::get().entropy.entropy_start_frame_num == i)
       {
         if (Config::get().entropy.entropy_use_internal)
         {
           matrix_aligned = transformToOneline(coords, Config::get().entropy.entropy_internal_dih, true);
         }
         else
         {
           matrix_aligned = transformToOneline(coords, Config::get().entropy.entropy_trunc_atoms_num, false);
         }
       }
       //Building one huge [coordinates] x [frames] matrix by appending for every frame

     }
   }
   transpose(matrix_aligned);
   //NECESSARY because of implementation details, don't worry about it for now; rows are DOFs, columns are frames FROM HERE ON!

   if (!Config::get().entropy.entropy_use_internal)
   {
     massweight(matrix_aligned, coords_ref, true);
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