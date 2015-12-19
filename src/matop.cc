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
    Matrix_Class out_mat((unsigned int)(in.size()), 3u);
    for (unsigned int l = 0; l < (unsigned int) in.size(); l++)
    {
      coords::cartesian_type tempcoord2;
      tempcoord2 = in.xyz(l);
      out_mat(l, 0) = tempcoord2.x();
      out_mat(l, 1) = tempcoord2.y();
      out_mat(l, 2) = tempcoord2.z();
    }
    return out_mat.transposed();
  };

  Matrix_Class transfer_to_matr_internal(coords::Coordinates const& in)
  {
    const unsigned int sizer = (unsigned int) in.size();
    Matrix_Class out_mat(sizer, 3u);
    for (size_t l = 0; l < in.size(); l++)
    {
      out_mat(l, 0) = in.intern(l).radius();
      out_mat(l, 1) = in.intern(l).inclination().radians();
      out_mat(l, 2) = in.intern(l).azimuth().radians();
    }
    return out_mat.transposed();
  };

  Matrix_Class transform_coordinates(coords::Coordinates& input)
  {
	  Matrix_Class output(3, (unsigned int) input.size());
	  for (size_t l = 0; l < (unsigned int) input.size(); l++)
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

    for (unsigned int i = 0; i < input.cols(); i++)
    {
      coords::Cartesian_Point tempcoord2(input(0, i), input(1, i), input(2, i));
      tempcoord1.push_back(tempcoord2);
    }
    return tempcoord1;
  }

  coords::Representation_Internal transfer_to_internalRepressentation(Matrix_Class const& input)
  {
    coords::Representation_Internal tempcoord1;

    for (unsigned int i = 0; i < input.cols(); i++)
    {
      coords::internal_type tempcoord2(input(0, i), scon::ang<float_type>(input(1, i)), scon::ang<float_type>(input(2, i)));
      tempcoord1.push_back(tempcoord2);
    }
    return tempcoord1;
  }

  Matrix_Class transform_3n_nf(Matrix_Class const& input)
  {
    Matrix_Class transformed_matrix(1, (input.rows() * input.cols()));
    int j = 0;
    for (unsigned int i = 0; i < input.cols(); i++)
    {
      for (unsigned int k = 0; k < input.rows(); k++)
      {
        transformed_matrix(0, j + k) = input(k, i);
      }
      j = j + input.rows();
    }
    return transformed_matrix;
  }

  void massweight(Matrix_Class& input, coords::Coordinates const& coords, bool to_meter, std::vector<unsigned int> atomsThatAreUsed)
  //coords are reference coord object to get atomic masses from
  //boolean controls whether coords should also be multiplied with 10e-10 to convert angstrom to meters)
  {
    if (atomsThatAreUsed.empty())
    {
      for (unsigned int i = 0; i < input.rows(); i = i + 3)
      {
        double temp = sqrt(coords.atoms(i / 3u).mass() * 1.6605402 * 10e-27);
        if (to_meter)
        {
          temp *= 10e-10;
        }
        for (unsigned int j = 0; j < input.cols(); j++)
        {
          for (unsigned int k = 0; k < 3; k++)
          {
            (input)(i + k, j) *= temp;
          }
        }
      }
    }
    else
    {
      for (unsigned int i = 0; i < input.rows(); i = i + 3)
      {
        double temp = sqrt(coords.atoms(atomsThatAreUsed[i / 3u]).mass() * 1.6605402 * 10e-27);
        if (to_meter)
        {
          temp *= 10e-10;
        }
        for (unsigned int j = 0; j < input.cols(); j++)
        {
          for (unsigned int k = 0; k < 3; k++)
          {
            (input)(i + k, j) *= temp;
          }
        }
      }
    }
  }

  void undoMassweight(Matrix_Class& input, coords::Coordinates const& coords, bool to_meter, std::vector<unsigned int> atomsThatAreUsed)
  {
    if (atomsThatAreUsed.empty())
    {
      for (unsigned int i = 0; i < input.rows(); i = i + 3)
      {
        double temp = sqrt(coords.atoms(i / 3u).mass() * 1.6605402 * 10e-27);
        if (to_meter)
        {
          temp *= 10e-10;
        }
        for (unsigned int j = 0; j < input.cols(); j++)
        {
          for (unsigned int k = 0; k < 3; k++)
          {
            (input)(i + k, j) /= temp;
          }
        }
      }
    }
    else
    {
      for (unsigned int i = 0; i < input.rows(); i = i + 3)
      {
        double temp = sqrt(coords.atoms(atomsThatAreUsed[i / 3u]).mass() * 1.6605402 * 10e-27);
        if (to_meter)
        {
          temp *= 10e-10;
        }
        for (unsigned int j = 0; j < input.cols(); j++)
        {
          for (unsigned int k = 0; k < 3; k++)
          {
            (input)(i + k, j) /= temp;
          }
        }
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
    void prepare_pca(Matrix_Class const& input, Matrix_Class& eigenvalues, Matrix_Class& eigenvectors, Matrix_Class& pca_modes, int rank)
    {
      Matrix_Class cov_matr = (input.transposed());
      Matrix_Class ones(input.cols(), input.cols());
      ones.fillwith(1.0);
      cov_matr = cov_matr - ones * cov_matr / input.cols();
      cov_matr = cov_matr.transposed() * cov_matr;
      cov_matr = cov_matr / input.cols();
      float_type cov_determ = 0.;
      int *cov_rank = new int;
      cov_matr.eigensym(eigenvalues, eigenvectors, cov_rank);
      auto debugg = cov_matr.to_std_vector();
      rank = *cov_rank;
      if (rank < (int)eigenvalues.rows() || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90))
      {
        std::cout << "Notice: covariance matrix is singular.\n";
        std::cout << "Details: rank of covariance matrix is " << *cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << ".\n";
        if (Config::get().PCA.pca_remove_dof)
        {
          int temp = std::max(6, int((cov_matr.rows() - *cov_rank)));
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
      for (unsigned int i = 0; i < eigenvalues.rows(); i++)
      {
        sum_of_all_variances += eigenvalues(i);
      }

      if (Config::get().PCA.pca_trunc_var != 1 && Config::get().PCA.pca_trunc_var < 1 && Config::get().PCA.pca_trunc_var > 0)
      {
        double temporary_sum_of_variances = 0.0;
        Matrix_Class pca_modes2 = pca_modes;
        if (eigenvalues.rows() < 2u) throw("Eigenvalues not initialized, error in truncating PCA values. You are in trouble.");
        for (unsigned int i = 0u - 1; i < eigenvalues.rows(); i++)
        {
          if (Config::get().PCA.pca_trunc_var < temporary_sum_of_variances / sum_of_all_variances)
          {
            int rows_of_pca_modes2 = pca_modes2.rows();
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

    Matrix_Class transform_3n_nf_trunc_pca(Matrix_Class const& input)
    {
      //Transform it as necessary
      Matrix_Class transformed_matrix(1u, (unsigned int)( (3 * Config::get().PCA.pca_trunc_atoms_num.size()) ));
      int j = 0;
      unsigned int quicksearch = 0;
      for (unsigned int i = 0; i < input.cols(); i++) //iterate over atoms
      {
        bool checker = false;
        for (unsigned int l = quicksearch; l < Config::get().PCA.pca_trunc_atoms_num.size(); l++) //iterate over vector of atoms to account for
        {
          if (Config::get().PCA.pca_trunc_atoms_num[l] - 1 == i)
          {
            checker = true;
            quicksearch++;
            break;
          }
        }
        if (checker)
        {
          for (unsigned int k = 0; k < input.rows(); k++)
          {
            transformed_matrix(0, j + k) = (input)(k, i);
          }
          j = j + input.rows();
        }
      }
      return transformed_matrix;
    }

    Matrix_Class transform_3n_nf_internal_pca(Matrix_Class const& input)
    {
      // HAS TO BE TESTED NOT TESTED YET but should be functional
      Matrix_Class out;
      Matrix_Class transformed_matrix(1u, (unsigned int)((Config::get().PCA.pca_internal_dih.size() * 2 + 2 * Config::get().PCA.pca_internal_ang.size() - Config::get().PCA.pca_internal_bnd.size())));
      int j = 0;
      unsigned int quicksearch_dih = 0;
      unsigned int quicksearch_ang = 0;
      unsigned int quicksearch_bnd = 0;
      for (unsigned int i = 0; i < input.cols(); i++)
      {
        unsigned int keeper = 0;
        bool checker_dih = false;
        bool checker_ang = false;
        bool checker_bnd = false;
        for (unsigned int l = quicksearch_dih; l < Config::get().PCA.pca_internal_dih.size(); l++)
        {
          if (Config::get().PCA.pca_internal_dih.size() != 0)
          {
            if (Config::get().PCA.pca_internal_dih[l] - 1 == i)
            {
              checker_dih = true;
              quicksearch_dih++;
              break;
            }
          }
        }
        for (unsigned int l = quicksearch_ang; l < Config::get().entropy.entropy_internal_ang.size(); l++)
        {
          if (Config::get().PCA.pca_internal_ang.size() != 0)
          {
            if (Config::get().PCA.pca_internal_ang[l] - 1 == i)
            {
              checker_ang = true;
              quicksearch_ang++;
              break;
            }
          }
        }
        for (unsigned int l = quicksearch_bnd; l < Config::get().entropy.entropy_internal_bnd.size(); l++)
        {
          if (Config::get().PCA.pca_internal_bnd.size() != 0)
          {
            if (Config::get().PCA.pca_internal_bnd[l] - 1 == i)
            {
              checker_bnd = true;
              quicksearch_bnd++;
              break;
            }
          }
        }
        if (checker_bnd)
        {
          transformed_matrix(0, j + keeper) = (input)(0, i);
          keeper++;
        }
        if (checker_ang)
        {
          transformed_matrix(0, j + keeper) = cos((input)(1, i));
          transformed_matrix(0, j + keeper + 1) = sin((input)(1, i));
          keeper += 2;
        }
        if (checker_dih)
        {
          transformed_matrix(0, j + keeper) = cos((input)(2, i));
          transformed_matrix(0, j + keeper + 1) = sin((input)(2, i));
          keeper += 2;
        }
        j = j + keeper;
      }
      out = transformed_matrix;
      return out;
    }

    void output_probability_density(Matrix_Class& pca_modes)
    {
      using namespace histo;
      std::vector<unsigned int> dimensionsToBeUsedInHistogramming = Config::get().PCA.pca_dimensions_for_histogramming;

      if (dimensionsToBeUsedInHistogramming.size() > 2)
      {
        std::cerr << "More than 2 dimensions for histogramming is not yet supported!\n";
        std::cerr << "Working with 2 dimensions";
        dimensionsToBeUsedInHistogramming.resize(2);
      }
      else if (dimensionsToBeUsedInHistogramming.size() == 1)
      {
        //DUMMY
        Histograms<float_type> *histograms_p;
        //Initializing
        if (Config::get().PCA.pca_histogram_number_of_bins > 0u)
        {
          size_t histogramBins = Config::get().PCA.pca_histogram_number_of_bins;
          histograms_p = new Histograms<float_type>((size_t) pca_modes.rows(), histogramBins);
        }
        else if (Config::get().PCA.pca_histogram_width > 0.)
        {
          float_type histogramWidth = Config::get().PCA.pca_histogram_number_of_bins;
          histograms_p = new Histograms<float_type>((size_t) pca_modes.rows(), histogramWidth);
        }
        else
        {
          throw "Error in output_probability_density, exiting.\n You need to specify either a numbger of histogram bins or a bin width in the inputfile.\n";
        }

        std::cout << "Starting: Histogramming...";
        //Filling the histogram
        for (unsigned int i = 0u; i < pca_modes.rows(); i++)
        {
          for (unsigned int j = 0u; j < pca_modes.cols(); j++)
          {
            histograms_p->add_value(i, pca_modes(i,j));
          }
        }
        histograms_p->distribute();
        std::cout << *histograms_p;
        delete histograms_p;
      }
      else if (dimensionsToBeUsedInHistogramming.size() == 2)
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
          float_type histogramWidth = Config::get().PCA.pca_histogram_number_of_bins;
          histograms_p = new DimensionalHistogram<float_type>( (size_t)dimensionsToBeUsedInHistogramming.size(), histogramWidth);
        }
        else
        {
          throw "Error in output_probability_density, exiting.\n You need to specify either a numbger of histogram bins or a bin width in the inputfile.\n";
        }

        std::cout << "Starting: Histogramming...";
        //Filling the histogram
        for (unsigned int j = 0u; j < pca_modes.cols(); j++)
        {
          std::vector<float_type> temp(dimensionsToBeUsedInHistogramming.size());
          for (unsigned int i = 0u; i < dimensionsToBeUsedInHistogramming.size(); i++)
          {
            temp[i] = pca_modes(dimensionsToBeUsedInHistogramming[i] - 1, j);
          }
          histograms_p->add_value(temp);
        }
        
        histograms_p->distribute();
        histograms_p->write("pca_histograms");
        histograms_p->writeAuxilaryData("pca_histograms_auxdata");
        delete histograms_p;
      }
    }

    void readEigenvectorsAndModes(Matrix_Class& eigenvectors, Matrix_Class& trajectory, std::string& additionalInformation, std::string filename)
    {
      std::ifstream pca_modes_stream(filename, std::ios::in);
      std::string line;
      std::getline(pca_modes_stream, line);
      int dimensions = std::stoi(line.substr(13, 2));
      while (line.find("Eigenvectors") == std::string::npos)
      {
        std::getline(pca_modes_stream, line);
      }

      std::getline(pca_modes_stream, line);
      eigenvectors.resize(stoi(line.substr(7, 10)), stoi(line.substr(20, 10)));

      for (unsigned int i = 0u; i < eigenvectors.rows(); i++)
      {
        std::getline(pca_modes_stream, line);
        size_t whitespace = 0u, lastWhitespace = 0u;
        for (unsigned int j = 0u; j < eigenvectors.cols(); j++)
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

      for (unsigned int i = 0u; i < trajectory.rows(); i++)
      {
        std::getline(pca_modes_stream, line);
        size_t whitespace = 0u, lastWhitespace = 0u;
        for (unsigned int j = 0u; j < trajectory.cols(); j++)
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

    float_type knn_distance(Matrix_Class const& input, unsigned int const& dimension_in, unsigned int const& k_in, unsigned int const& row_queryPt, unsigned int const& col_queryPt)
    {
      float_type temp_distance = 0.0;
      float_type hold_distance;
      float_type *distanceList = new float_type[k_in];
      distanceList[0] = std::numeric_limits<float_type>::max();

      // This iterates over the "n-th" next neighbors
      // (to get the second next neighbor you have to find the first next neighbor etc. )
      for (unsigned int i = 0; i < k_in; i++)
      {
        // Get max() is initial value for distance comparison.
        // This garantues that the first calcualted value is smaller than
        // hold_distance.
        hold_distance = std::numeric_limits<float_type>::max();

        // This iterates over all points in the set
        for (unsigned int j = 0; j < input.cols(); j++)
        {
          // Of course we cannot count the distance of an member to itself
          // since it is =0.0
          if (j == col_queryPt) { continue; }

          // For number of dimensions, add the squared distance of the queryPt
          // to the current point ("j") of each dimensions which equals a 
          // squared distance in euclidean space
          temp_distance = 0.0;
          for (unsigned int l = 0; l < dimension_in; l++)
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

    float_type knn_distance(Matrix_Class const& input, unsigned int const& dimension_in, unsigned int const& k_in, std::vector<unsigned int>& row_queryPts, unsigned int const& col_queryPt)
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
      for (unsigned int i = 0; i < k_in; i++)
      {
        // Get max() is initial value for distance comparison.
        // This garantues that the first calcualted value is smaller than
        // hold_distance.
        hold_distance = std::numeric_limits<float_type>::max();

        for (unsigned int j = 0; j < input.cols(); j++)
        {
          // Of course we cannot count the distance of an member to itself
          // since it is =0.0
          if (j == col_queryPt) { continue; }

          // For number of dimensions, add the squared distance of the queryPt
          // to the current point ("j") of each dimensions which equals a 
          // squared distance in euclidean space
          temp_distance = 0.0;
          for (unsigned int l = 0; l < dimension_in; l++)
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
      Matrix_Class cov_matr = (input.transposed());
      cov_matr = cov_matr - Matrix_Class( input.cols(), input.cols() ).filledwith(1.) * cov_matr / input.cols();
      cov_matr = cov_matr.transposed() * cov_matr;
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
          unsigned int temp = (unsigned int)std::max(6, int((cov_matr.rows() - *cov_rank)));
		      eigenvalues.shed_rows((unsigned int) eigenvalues.rows() - temp, (unsigned int) eigenvalues.rows() - 1u);
		      eigenvectors.shed_cols((unsigned int) eigenvectors.cols() - temp, (unsigned int) eigenvectors.cols() - 1u);
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
      for (unsigned int i = 0; i < int(eigenvalues.rows()); i++)
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
      Matrix_Class eigenvectors_t(eigenvectors.transposed());
      Matrix_Class pca_modes = eigenvectors_t * input;
      Matrix_Class entropy_anharmonic(pca_modes.rows());
      entropy_anharmonic.fillwith(0.0);
      Matrix_Class entropy_mi(pca_modes.rows(), pca_modes.rows());
      entropy_mi.fillwith(0.0);
      Matrix_Class classical_entropy(pca_modes.rows());
      classical_entropy.fillwith(0.0);
      int const size = entropy_anharmonic.rows();


      // II. Calculate Non-Paramteric Entropies
      // Marginal
      for (unsigned int i = 0; i < pca_modes.rows(); i++)
      {
        float_type distance = 0.0;
        for (unsigned int k = 0; k < pca_modes.cols(); k++)
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
          for (unsigned int k = 1; k < Config::get().entropy.entropy_method_knn_k; k++)
          {
            temp += 1.0 / float_type(i);
          }
        }
        distance -= temp;
        distance += 0.5772156649015328606065;
        entropy_anharmonic(i) = distance;
        distance = 0;

        //MI
        for (unsigned int j = i + 1; j < pca_modes.rows(); j++)
        {
          std::vector<unsigned int> query_rows{ i,j };
          for (unsigned int k = 0; k < pca_modes.cols(); k++)
          {
            distance += log(sqrt(knn_distance(pca_modes, 2, Config::get().entropy.entropy_method_knn_k, query_rows, k)));
          }

          distance /= 2 * float_type(pca_modes.cols());
          float_type temp2 = float_type(pca_modes.cols()) * pow(3.14159265358979323846, 1);

          distance += log(temp2);
          temp2 = 0;
          if (Config::get().entropy.entropy_method_knn_k != 1)
          {
            for (unsigned int u = 1; u < Config::get().entropy.entropy_method_knn_k; u++)
            {
              temp2 += 1.0 / float_type(u);
            }
          }
          distance -= temp2;
          distance += 0.5772156649015328606065;
          entropy_mi(i, j) = distance;

        }
      }

      unsigned int counterForLargeNegativeM_I_Terms = 0u;
      for (unsigned int i = 0; i < entropy_anharmonic.rows(); i++)
      {
        for (unsigned int j = (i + 1); j < entropy_anharmonic.rows(); j++)
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
      for (unsigned int i = 0; i < entropy_anharmonic.rows(); i++)
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
      for (unsigned int i = 0; i < entropy_anharmonic.rows(); i++)
      {
        delta_entropy += entropy_anharmonic(i);
        for (unsigned int j = (i + 1); j < entropy_anharmonic.rows(); j++)
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
      Matrix_Class marginal_entropy_storage(input.rows(), 1u);
      marginal_entropy_storage.fillwith(0.0);

      float_type distance = knn_distance(input, input.rows(), Config::get().entropy.entropy_method_knn_k, 0, 1);
      for (unsigned int k = 1; k < input.cols(); k++)
      {
        //Should be equivalent to the original version but with less costly sqrt:
        distance *= knn_distance(input, input.rows(), Config::get().entropy.entropy_method_knn_k, 0, k);
      }
      distance = std::log(std::sqrt(distance));

      // Calculate Non-Paramteric Entropies (old)
      /*float_type distance = 0.0;
      for (unsigned int k = 0; k < this->cols(); k++)
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
        for (unsigned int i = 1; i < Config::get().entropy.entropy_method_knn_k; i++)
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
      Matrix_Class marginal_entropy_storage(input.rows());
      marginal_entropy_storage.fillwith(0.0);

      //Calculate Non-Paramteric Entropies
      for (unsigned int i = 0; i < input.rows(); i++)
      {
        double distance = 0.0;
        for (unsigned int k = 0; k < input.cols(); k++)
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
          for (unsigned int k = 1; k < Config::get().entropy.entropy_method_knn_k; k++)
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
      for (unsigned int i = 0; i < marginal_entropy_storage.rows(); i++)
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
      Matrix_Class cov_matr = (input.transposed());
      cov_matr = cov_matr - Matrix_Class(input.cols(), input.cols()).filledwith(1.) * cov_matr / (float_type)input.cols();
      cov_matr = cov_matr.transposed() * cov_matr;
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
          unsigned int temp = (unsigned int)std::max(6, int((cov_matr.rows() - *cov_rank)));
          eigenvalues.shed_rows((unsigned int)(eigenvalues.rows()) - temp, eigenvalues.rows() - 1u);
          eigenvectors.shed_cols((unsigned int)(eigenvectors.cols()) - temp, eigenvectors.cols() - 1u);
        }
        else
        {
		      eigenvalues.shed_rows( (*cov_rank), (unsigned int) (eigenvalues.rows()) - 1u);
	    	  eigenvectors.shed_cols( (*cov_rank), (unsigned int) (eigenvectors.cols()) - 1u);
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
      for (unsigned int i = 0; i < eigenvalues.rows(); i++)
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
      Matrix_Class cov_matr = (input.transposed());
      cov_matr = cov_matr - Matrix_Class( input.cols(), input.rows() ).filledwith(1.0) * cov_matr / input.cols();
      cov_matr = cov_matr.transposed() * cov_matr;
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
      Matrix_Class cov_matr = (input.transposed());
      cov_matr = cov_matr - Matrix_Class(input.cols(), input.cols()).filledwith(1.) * cov_matr / input.cols();
      cov_matr = cov_matr.transposed() * cov_matr;
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
      cov_matr = cov_matr + Matrix_Class(cov_matr.rows(), cov_matr.cols()).return_identity();
      float_type entropy_sho = cov_matr.determ();

      entropy_sho = log(entropy_sho) * 0.5 * 1.38064813 * 6.02214129 * 0.239;
      //This stems from: k_B * Avogadro * (to_calories) *0.5

      std::cout << "Entropy in QH-approximation: " << entropy_sho << " cal / (mol * K)" << std::endl;
      return entropy_sho;
    }
	
    Matrix_Class transform_3n_nf_trunc_entropy(Matrix_Class const& input)
    { // Transform it as necessary
      Matrix_Class transformed_matrix(1u, (unsigned int)((input.rows() * input.cols() - 3 * Config::get().entropy.entropy_trunc_atoms_num.size())));
      int j = 0;
      unsigned int quicksearch = 0;
      for (unsigned int i = 0; i < input.cols(); i++)
      {
        bool checker = false;
        for (unsigned int l = quicksearch; l < Config::get().PCA.pca_trunc_atoms_num.size(); l++)
        {
          if (Config::get().PCA.pca_trunc_atoms_num[l] - 1 == i)
          {
            checker = true;
            quicksearch++;
            break;
          }
        }
        if (checker)
        {
          for (unsigned int k = 0; k < input.rows(); k++)
          {
            transformed_matrix(0, j + k) = (input)(k, i);
          }
        }
        j = j + input.rows();
      }
      return transformed_matrix;
    }

    Matrix_Class transform_3n_nf_internal_entropy(Matrix_Class const& input)
    {
      // HAS TO BE TESTED NOT TESTED YET but should be functional
      Matrix_Class out;
      Matrix_Class transformed_matrix(1u, (unsigned int)((Config::get().entropy.entropy_internal_dih.size() * 2 + Config::get().entropy.entropy_internal_ang.size() * 2 + Config::get().entropy.entropy_internal_bnd.size())));
      int j = 0;
      unsigned int quicksearch_dih = 0;
      unsigned int quicksearch_ang = 0;
      unsigned int quicksearch_bnd = 0;
      for (unsigned int i = 0; i < input.cols(); i++)
      {
        unsigned int keeper = 0;
        bool checker_dih = false;
        bool checker_ang = false;
        bool checker_bnd = false;
        if (Config::get().entropy.entropy_internal_dih.size() != 0)
        {
          for (unsigned int l = quicksearch_dih; l < Config::get().entropy.entropy_internal_dih.size(); l++)
          {

            if (Config::get().entropy.entropy_internal_dih[l] - 1 == i)
            {
              checker_dih = true;
              quicksearch_dih++;
              break;
            }
          }
        }
        for (unsigned int l = quicksearch_ang; l < Config::get().entropy.entropy_internal_ang.size(); l++)
        {
          if (Config::get().entropy.entropy_internal_ang.size() != 0)
          {
            if (Config::get().entropy.entropy_internal_ang[l] - 1 == i)
            {
              checker_ang = true;
              quicksearch_ang++;
              break;
            }
          }
        }
        for (unsigned int l = quicksearch_bnd; l < Config::get().entropy.entropy_internal_bnd.size(); l++)
        {
          if (Config::get().entropy.entropy_internal_bnd.size() != 0)
          {
            if (Config::get().entropy.entropy_internal_bnd[l] - 1 == i)
            {
              checker_bnd = true;
              quicksearch_bnd++;
              break;
            }
          }
        }
        if (checker_bnd)
        {
          transformed_matrix(0, j + keeper) = (input)(0, i);
          keeper++;
        }
        if (checker_ang)
        {
          transformed_matrix(0, j + keeper) = cos((input)(1, i));
          transformed_matrix(0, j + keeper + 1) = sin((input)(1, i));
          keeper += 2;
        }
        if (checker_dih)
        {
          transformed_matrix(0, j + keeper) = cos((input)(2, i));
          transformed_matrix(0, j + keeper + 1) = sin((input)(2, i));
          keeper += 2;
        }
        j = j + keeper;
      }
      out = transformed_matrix;
      return out;
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
    float_type calc_rmsd(Matrix_Class const& input, Matrix_Class const& ref)
      //LAYER 1
    {
      unsigned int size = input.cols();
      float_type temp(0.), temp2(0.), temp3(0.);
      for (unsigned int i = 0; i < size; i++) {
        temp = 0;
        temp3 = 0;
        for (unsigned int j = 0; j < input.cols(); j++) {
          temp = (ref(j, i) - (input)(j, i));
          temp *= temp;
          temp3 += temp;
        }
        temp2 += temp3;
      }
      return (sqrt(temp2 / float_type(size)));
    }

    float_type drmsd_calc(Matrix_Class const& input, Matrix_Class const& ref)
      //LAYER 1
    {
      float_type value = 0;
      int counter(0);
      Matrix_Class matr_structure(input);
      for (unsigned int i = 0; i < input.cols(); i++) {
        for (unsigned int j = 0; j < i; j++)
        {
          float_type holder = sqrt((ref(0, i) - ref(0, j)) * (ref(0, i) - ref(0, j)) + (ref(1, i) - ref(1, j))* (ref(1, i) - ref(1, j)) + (ref(2, i) - ref(2, j))*(ref(2, i) - ref(2, j)));
          float_type holder2 = sqrt(((input)(0, i) - (input)(0, j)) * ((input)(0, i) - (input)(0, j)) + ((input)(1, i) - (input)(1, j)) * ((input)(1, i) - (input)(1, j)) + ((input)(2, i) - (input)(2, j)) * ((input)(2, i) - (input)(2, j)));
          value += (holder2 - holder) * (holder2 - holder);
          counter++;
        }
      }
      return sqrt(value / counter);
    }

    float_type holmsander_calc(Matrix_Class const& input, Matrix_Class const& ref)
      //LAYER 1
    {
      float_type value(0);
      Matrix_Class matr_structure = input;
      for (unsigned int i = 0; i < input.cols(); i++) {
        for (unsigned int j = 0; j < i; j++)
        {
          float_type holder = sqrt((ref(0, i) - ref(0, j)) * (ref(0, i) - ref(0, j)) + (ref(1, i) - ref(1, j))* (ref(1, i) - ref(1, j)) + (ref(2, i) - ref(2, j))*(ref(2, i) - ref(2, j)));
          float_type holder2 = sqrt(((input)(0, i) - (input)(0, j)) * ((input)(0, i) - (input)(0, j)) + ((input)(1, i) - (input)(1, j))* ((input)(1, i) - (input)(1, j)) + ((input)(2, i) - (input)(2, j))*((input)(2, i) - (input)(2, j)));
          value += abs(holder2 - holder) * exp(-1 * (holder2 + holder)*(holder2 + holder) / (4 * (config::align().holm_sand_r0) * (config::align().holm_sand_r0))) / (holder2 + holder);
        }
      }
      return value;
    }

    Matrix_Class rotated(Matrix_Class const& input, Matrix_Class const& ref)
      //LAYER 1
    {
      Matrix_Class working_copy(input, false);
      rotate(working_copy, ref);
      return working_copy;
    }

    void rotate(Matrix_Class& input, Matrix_Class const& ref)
      //LAYER 1
    {
      Matrix_Class c(input * ref.transposed());
      //Creates Covariance Matrix

      Matrix_Class s, V, U;
      c.singular_value_decomposition(U, s, V);

      Matrix_Class unit(c.rows(), c.rows()); // Create empty dummy matrix of right size for unitary
      unit.identity(); // Make identity matrix
      if ((c.det_sign() < 0)) //Making sure that U will do a proper rotation (rows/columns have to be right handed system)
      {
        unit(2, 2) = -1;
      }
      U.transpose();
      unit = unit * U;
      unit = V * unit;
      input = unit * input;
    }

    void align_center_of_mass(Matrix_Class& input, coords::Coordinates const& coords_in)
      //LAYER 1
    {
      coords::Cartesian_Point com_ref = coords_in.center_of_mass();

      for (unsigned int i = 0; i < input.cols(); i++)
      {
        (input)(0, i) -= com_ref.x();
        (input)(1, i) -= com_ref.y();
        (input)(2, i) -= com_ref.z();
      }
    }

    void align_center_of_geo(Matrix_Class& input)
      //LAYER 1
    {
      std::cerr << "CAUTION!!!!!! A DEPRECATED FUNCTION IS USED! you used Matrix_Class::align(void) which alignes to center of geometry, better use Matrix_Class::align(coords::coordinates&) which alignes to center of mass! This is the usual procedure!\n";
      std::vector <float_type> com_ref = center_of_geo(input);
      for (unsigned int j = 0; j < input.rows(); j++) {
        for (unsigned int i = 0; i < input.cols(); i++) {
          (input)(j, i) -= com_ref[j];
        }
      }
    }

    std::vector <float_type> center_of_geo(Matrix_Class const& input)
      //LAYER 1
    {
      std::vector <float_type> output(3, 0.0);
      float_type mean = 0;
      for (unsigned int j = 0; j < 3; j++) {
        mean = 0.0;
        for (unsigned int i = 0; i < input.rows(); i++) {
          mean += (input)(i, j);
        }
        output[j] = (mean / (float_type(input.rows())));
      }
      return output;
    }
  }

} //END NAMESPACE matop