#include "entropy.h"
#pragma once


namespace entropy
{

  float_type knn_distance_eucl_squared(Matrix_Class const& input, size_t const& dimension_in, size_t const& k_in, std::vector<size_t> const& row_queryPts, size_t const& col_queryPt, coords::float_type* buffer)
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


  float_type maximum_norm_knn_distance(Matrix_Class const& input, size_t const& dimension_in, size_t const& k_in, std::vector<size_t> const& row_queryPts, size_t const& col_queryPt, coords::float_type* buffer)
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

        // Maximum norm evaluation
        temp_distance = 0.0;
        for (size_t l = 0; l < dimension_in; l++)
        {
          temp_distance = std::max(temp_distance, std::abs((input)(row_queryPts[l], j) - (input)(row_queryPts[l], col_queryPt)));
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

  float_type maximum_norm_knn_distance_with_pointiters(Matrix_Class const& input, size_t const& dimension_in, size_t const& k_in, std::vector<size_t> const& row_queryPts,
    size_t const& col_queryPt, coords::float_type* buffer, size_t* intbuffer)
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
    size_t *intbfr = intbuffer;
    if (buffer == nullptr)
    {
      distanceList = new float_type[k_in];
    }
    if (intbuffer == nullptr)
    {
      throw std::runtime_error("intbuffer cannot be nullptr for lombardi entropy procedure.");
    }
    intbfr[0] = std::numeric_limits<size_t>::quiet_NaN();
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

        // Maximum norm evaluation
        temp_distance = 0.0;
        for (size_t l = 0; l < dimension_in; l++)
        {
          temp_distance = std::max(temp_distance, std::abs((input)(row_queryPts[l], j) - (input)(row_queryPts[l], col_queryPt)));
        }

        // If we are searching for the actual, "first" nearest neighbor
        // we will compare the distance in "hold_distance", initialized as a huge number,
        // to the calculated temp_distance ("Is this point nearer to the query point
        // than the previously established nearest point?"). We will, of course,
        // keep the smaller (squared) distance of the two.
        if (i == 0 && temp_distance < hold_distance)
        {
          hold_distance = temp_distance;
          intbuffer[i] = j;
        }

        // If we are searching for the "(i + 1)-th" nearest neighbot
        // we will compare "hold_distance" to "temp_distance". If temp_distance is smaller, it
        // will be kept, however, only if it is also larger than the previously established
        // "i-th" nearest neighbor. Otherwise we would always just get the absolute nearest neighbor in
        // a set of points, and not the k-th nearest neighbor.
        else if (i > 0 && temp_distance < hold_distance && temp_distance > distanceList[i - 1])
        {
          hold_distance = temp_distance;
          intbuffer[i] = j;
        }
      }
      distanceList[i] = hold_distance;
    }
    float_type keeper = distanceList[k_in - 1u];
    if (buffer == nullptr) delete[] distanceList;
    return keeper;
  }


  TrajectoryMatrixRepresentation::TrajectoryMatrixRepresentation(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords)
  {
    generateCoordinateMatrix(ci, coords);
  };


  void TrajectoryMatrixRepresentation::generateCoordinateMatrix(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords)
  {
    // First, adjust number of truncated atoms to be used to zero, 
    // in case truncation is not to be used (duh)
    if (!Config::get().entropy.entropy_trunc_atoms_bool) 
      Config::set().entropy.entropy_trunc_atoms_num = std::vector<size_t>();

    // Initialize the reference frame (for alignment etc)
    coords::Coordinates coords_ref(coords);
    auto holder = (*ci).PES()[Config::get().entropy.entropy_ref_frame_num].structure.cartesian;
    coords_ref.set_xyz(holder);



    //////////////
    //
    //  A L I G N M E N T
    //
    if (Config::get().entropy.entropy_alignment && Config::get().entropy.entropy_use_internal)
    {
      std::cout << "Alignment is (in this case) redundant since internal coordinates are used. Alignment is skipped. Check your INPUTFILE please.\n";
      std::cout << "Continuing anyway...\n";
    }
    // Translational alignment of the reference frame
    if (Config::get().entropy.entropy_alignment)
    {
      align::centerOfMassAlignment(coords_ref);
    }
    //
    //
    //////////////

    // The coordinateMatrix that will now
    // be built up step by step

    // Build up coordinate Matrix according to input file
    // This section will only set the right size of
    // the coordsMatrix object
    //
    // If internals are used a nonlinear transform is applied
    // according to Knapp (DOI 10.1063/1.2746330) to transform
    // the angular space to a linear space
    if (Config::get().entropy.entropy_use_internal)
    {
      coordsMatrix = Matrix_Class((size_t) /* explicitly casting to round down */ ((ci->size() - Config::get().entropy.entropy_start_frame_num) / Config::get().entropy.entropy_offset), \
        Config::get().entropy.entropy_internal_dih.size() * 2u);
    }
    // If truncated cartesians are desired, this section will
    // handle it.
    else if (Config::get().entropy.entropy_trunc_atoms_bool)
    {
      coordsMatrix = Matrix_Class((size_t) /* explicitly casting to round down */ ((ci->size() - Config::get().entropy.entropy_start_frame_num) / Config::get().entropy.entropy_offset), \
        Config::get().entropy.entropy_trunc_atoms_num.size() * 3u);
    }
    // This section is used if *all* cartesians of all atoms are used
    else
    {
      coordsMatrix = Matrix_Class((size_t) /* explicitly casting to round down */ ((ci->size() - Config::get().entropy.entropy_start_frame_num) / Config::get().entropy.entropy_offset), \
        coords.atoms().size() * 3u);
    }

    // Now the coordsMatrix will be filled with the coordinates
    // read in from ci
    //
    // j counts the (truncated) matrix access, i the frames in ci
    {
      size_t j = 0;
      // This section if internals are used
      // Remember: 
      // If internals are used a nonlinear transform is applied
      // according to Knapp (DOI 10.1063/1.2746330) to transform
      // the angular space to a linear space
      if (Config::get().entropy.entropy_use_internal)
      {
        for (size_t i = Config::get().entropy.entropy_start_frame_num; j < coordsMatrix.rows(); ++j, i += Config::get().entropy.entropy_offset)
        {
          auto holder2 = ci->PES()[i].structure.cartesian;
          coords.set_xyz(holder2);
          coordsMatrix.row(j) = ::matop::transformToOneline(coords, Config::get().entropy.entropy_internal_dih, true);
        }
      }
      // This section if cartesian coordinates are used
      // They will *later* be massweighted
      else
      {
        for (size_t i = Config::get().entropy.entropy_start_frame_num; j < coordsMatrix.rows(); ++j, i += Config::get().entropy.entropy_offset)
        {
          auto holder2 = ci->PES()[i].structure.cartesian;
          coords.set_xyz(holder2);
          // Translational and rotational alignment
          if (Config::get().entropy.entropy_alignment)
          {
            // Alignes center of mass
            align::centerOfMassAlignment(coords); 
            // Rotational alignment
            align::kabschAlignment(coords, coords_ref); 
          }
          coordsMatrix.row(j) = ::matop::transformToOneline(coords, Config::get().entropy.entropy_trunc_atoms_num, false);
        }
      }
    }

    // This transpose call is necessary
    // because of (legacy)implementation details, 
    // don't worry about it for now.
    // 
    // Rows are DOFs, Columns are frames FROM HERE ON!
    transpose(coordsMatrix);

    // Mass-weightening cartesian coordinates
    if (!Config::get().entropy.entropy_use_internal)
    {
      if (!Config::get().entropy.entropy_trunc_atoms_bool)
      {
        ::matop::massweight(coordsMatrix, coords_ref, true);
      }
      else
      {
        ::matop::massweight(coordsMatrix, coords_ref, true, Config::get().entropy.entropy_trunc_atoms_num);
      }
    }
  }

  float_type TrajectoryMatrixRepresentation::hnizdo(size_t const kNN)
  {
    // If openMP is used, every thread needs its own copy
    // Otherwise we may work with a reference
#ifdef _OPENMP
      Matrix_Class input(coordsMatrix);
#else
      Matrix_Class const& input(coordsMatrix);
#endif
      std::cout << "\nCommencing entropy calculation:\nNearest-Neighbor Nonparametric Method, according to Hnizdo et al. (DOI: 10.1002/jcc.20589)" << std::endl;
      
      std::vector<size_t> row_queryPoints;
      for (unsigned int i = 0u; i < input.rows(); ++i)
        row_queryPoints.push_back(i);
      


      float_type distance = 0.;
#ifdef _OPENMP
#pragma omp parallel firstprivate(input, kNN) reduction(+:distance)
      {
#endif
        KahanAccumulation<float_type> kahan_acc;
        float_type* buffer = new float_type[kNN];
#ifdef _OPENMP
#pragma omp for
#endif
        for (int i = 1; i < (int)input.cols(); i++)
        {
          // NOTE: This should be reviewed!
          //
          // Should be equivalent to the original version but with less costly sqrt:
          // This is the squared distance btw
          // distance *= knn_distance(input_workobj,  input.rows(), k, (size_t) 0u, (size_t) k);
          float_type distance = std::log(std::sqrt(knn_distance_eucl_squared(input, input.rows(), kNN, row_queryPoints, (size_t)i, buffer)));
        
          kahan_acc = KahanSum(kahan_acc, distance);
        }
        distance += kahan_acc.sum;
        delete[] buffer;
#ifdef _OPENMP
      }
#endif


      float_type hnizdoSum = distance;
      hnizdoSum /= float_type(input.cols()); //Number of Draws
      hnizdoSum *= float_type(input.rows()); //DIMENSIONS
      hnizdoSum += (log(float_type(input.cols())
        * pow(3.1415926535897932384626433832795029, float_type(float_type(input.rows())) / 2.) / (tgamma(0.5 * input.rows() + 1))));


      float_type temporary = 0;
      if (kNN != 1)
      {
        for (size_t i = 1; i < kNN; i++)
        {
          temporary += 1.0 / float_type(i);
        }
      }

      hnizdoSum -= temporary;
      hnizdoSum += 0.5772156649015328606065;

      float_type entropy = hnizdoSum;
      std::cout << "entropy (dimensionless): " << entropy << '\n';
      std::cout << "entropy: " << entropy * 1.380648813 * 6.02214129 * 0.239005736 << " cal / (mol * K)" << "\n";
      return entropy;
  }

  float_type TrajectoryMatrixRepresentation::hnizdo_marginal(size_t const kNN)
  {
    // If openMP is used, every thread needs its own copy
    // Otherwise we may work with a reference
#ifdef _OPENMP
    Matrix_Class input(coordsMatrix);
#else
    Matrix_Class const& input(coordsMatrix);
#endif // _OPENMP
    std::cout << "\nCommencing entropy calculation:\nNearest-Neighbor Nonparametric Method - only calculate sum of Marginal Entropies, according to Hnizdo et. al. (DOI: 10.1002/jcc.20589)" << std::endl;
    
    // Here the marginal entropies for each DOF are stored
    Matrix_Class marginal_entropy_storage(input.rows(), 1u, 0u);

    //Calculate Non-Paramteric Entropies
#ifdef _OPENMP
#pragma omp parallel for firstprivate(input, kNN) shared(marginal_entropy_storage)
#endif // _OPENMP
    for (int i = 0; i < (int)input.rows(); i++)
    {
      KahanAccumulation<double> kahan_acc;

      float_type* buffer = new float_type[kNN];
      for (size_t j = 0; j < input.cols(); j++)
      {
        float_type distance = log(sqrt(knn_distance_eucl_squared(input, 1, kNN, std::vector<size_t>{static_cast<size_t>(i)}, j, buffer))); // set eucledean distance to ouptut
        kahan_acc = KahanSum(kahan_acc, distance);
      }
      float_type distance = kahan_acc.sum / float_type(input.cols());

      // Adding mathematical constants according to original publication
      float_type temp = float_type(input.cols()) * pow(3.14159265358979323846, 0.5);
      temp /= tgamma((1. / 2.) + 1.);
      distance += log(temp);
      temp = 0;
      if (kNN != 1)
      {
        for (size_t l = 1; l < kNN; l++)
        {
          temp += 1.0 / float_type(l);
        }
      }
      distance -= temp;
      distance += 0.5772156649015328606065;
      marginal_entropy_storage(i) = distance;
      delete[] buffer;
    }

    //Calculate sum of of Entropies
    float_type entropy = 0;
    for (size_t i = 0; i < marginal_entropy_storage.rows(); i++)
    {
      entropy += marginal_entropy_storage(i);
    }
    std::cout << "entropy (dimensionless): " << entropy << '\n';
    std::cout << "entropy: " << entropy * 1.380648813 * 6.02214129 * 0.239005736 << " cal / (mol * K)" << "\n";
    return entropy;
  }

  float_type TrajectoryMatrixRepresentation::karplus()
  {
    std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Karplus et. al. (DOI 10.1021/ma50003a019)" << std::endl;
    Matrix_Class cov_matr = (transposed(coordsMatrix));
    cov_matr = cov_matr - Matrix_Class(coordsMatrix.cols(), coordsMatrix.cols(), 1.) * cov_matr / static_cast<float_type>(coordsMatrix.cols());
    cov_matr = transposed(cov_matr) * cov_matr;
    cov_matr = cov_matr / static_cast<float_type>(coordsMatrix.cols());
    float_type entropy = 0.0, cov_determ;
    if (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90)
    {
      std::cout << "Error: Covariance Matrix is singular. Try: a.) using internal coordinates b.) higher sampling rate in MD-Simulation c.) using more advanced methods.\n";
    }
    else if (cov_determ != 0)
    {
      entropy = log(pow(2 * SCON_PI, int(coordsMatrix.rows())) * cov_determ);
      entropy += (coordsMatrix.rows());
      entropy *= 1.380648813 * 6.02214129 * 0.239005736;
      std::cout << "Entropy in classical QH-approximation: " << entropy << " cal / (mol * K)" << std::endl;
    }
    return entropy;
  }

  float_type TrajectoryMatrixRepresentation::schlitter(float_type const temperatureInKelvin)
  {
    std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Schlitter (see: doi:10.1016/0009-2614(93)89366-P)" << std::endl;
    Matrix_Class cov_matr = transposed(coordsMatrix);
    cov_matr = cov_matr - Matrix_Class(coordsMatrix.cols(), coordsMatrix.cols(), 1.0) * cov_matr / static_cast<float_type>(coordsMatrix.cols());
    cov_matr = transposed(cov_matr) * cov_matr;
    cov_matr = cov_matr / static_cast<float_type>(coordsMatrix.cols());

    cov_matr *= (1.38064813 * /* 10e-23 J/K */ temperatureInKelvin * 2.718281828459 * 2.718281828459 / (1.054571726 /* * 10^-34 Js */ * 1.054571726 * 10e-45));
    cov_matr = cov_matr + Matrix_Class::Identity(cov_matr.rows(), cov_matr.cols());
    float_type entropy_sho = cov_matr.determ();

    entropy_sho = log(entropy_sho) * 0.5 * 1.38064813 * 6.02214129 * 0.239;
    //This stems from: k_B * Avogadro * (to_calories) *0.5

    std::cout << "Entropy in QH-approximation: " << entropy_sho << " cal / (mol * K)" << std::endl;
    return entropy_sho;
  }

  float_type TrajectoryMatrixRepresentation::knapp_marginal(float_type const temperatureInKelvin, bool removeDOF)
  {
    std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Knapp et. al. without corrections (Genome Inform. 2007;18:192-205.)" << std::endl;
    Matrix_Class cov_matr = (transposed(coordsMatrix));
    cov_matr = cov_matr - Matrix_Class(coordsMatrix.cols(), coordsMatrix.cols(), 1.) * cov_matr / (float_type)coordsMatrix.cols();
    cov_matr = transposed(cov_matr) * cov_matr;
    cov_matr *= (1.f / static_cast<float_type>(coordsMatrix.cols()));
    Matrix_Class eigenvalues;
    Matrix_Class eigenvectors;
    float_type cov_determ = 0.;
    int *cov_rank = new int;
    cov_matr.eigensym(eigenvalues, eigenvectors, cov_rank);

    //Remove Eigenvalues that should be zero if cov_matr is singular
    if ((*cov_rank < (int)eigenvalues.rows()) || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90))
    {
      std::cout << "Notice: covariance matrix is singular, attempting to fix by truncation of Eigenvalues.\n";
      std::cout << "Details: rank of covariance matrix is " << *cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << ".\n";
      if (removeDOF)
      {
        size_t temp = std::max(6, int((cov_matr.rows() - *cov_rank)));
        eigenvalues.shed_rows((eigenvalues.rows()) - temp, eigenvalues.rows() - 1u);
        eigenvectors.shed_cols((eigenvectors.cols()) - temp, eigenvectors.cols() - 1u);
      }
      else
      {
        eigenvalues.shed_rows((*cov_rank), (eigenvalues.rows()) - 1u);
        eigenvectors.shed_cols((*cov_rank), (eigenvectors.cols()) - 1u);
      }
    }
    else if (removeDOF)
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
    for (size_t i = 0; i < eigenvalues.rows(); i++)
    {
      pca_frequencies(i, 0u) = sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp / eigenvalues(i, 0u));
      alpha_i(i, 0u) = 1.05457172647 * 10e-34 / (sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp) * sqrt(eigenvalues(i, 0u)));
      quantum_entropy(i, 0u) = ((alpha_i(i, 0u) / (exp(alpha_i(i, 0u)) - 1)) - log(1 - exp(-1 * alpha_i(i, 0u)))) * 1.380648813 * 6.02214129 * 0.239005736;
      entropy_sho += quantum_entropy(i, 0u);
    }
    std::cout << "Entropy in QH-approximation: " << entropy_sho << " cal / (mol * K)" << std::endl;
    return entropy_sho;
  }

  float_type TrajectoryMatrixRepresentation::knapp(float_type const temperatureInKelvin, size_t const k_kNN, bool removeDOF)
  {
#ifdef _OPENMP
    Matrix_Class input(coordsMatrix);
#else
    Matrix_Class const& input(coordsMatrix);
#endif // _OPENMP

    std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Knapp et. al. with corrections (Genome Inform. 2007;18:192-205.)" << std::endl;
    Matrix_Class cov_matr = Matrix_Class{ transposed(input) };
    cov_matr = cov_matr - Matrix_Class(input.cols(), input.cols(), 1.) * cov_matr / static_cast<float_type>(input.cols());
    cov_matr = transposed(cov_matr) * cov_matr;
    cov_matr *= (1.f / static_cast<float_type>(input.cols()));
    Matrix_Class eigenvalues;
    Matrix_Class eigenvectors;
    float_type cov_determ = 0.;
    int *cov_rank = new int;
    cov_matr.eigensym(eigenvalues, eigenvectors, cov_rank);

    //Remove Eigenvalues that should be zero if cov_matr is singular
    if ((*cov_rank < (int)eigenvalues.rows()) || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90))
    {
      std::cout << "Notice: covariance matrix is singular, attempting to fix by truncation of Eigenvalues.\n";
      std::cout << "Details: rank of covariance matrix is " << *cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << ".\n";
      if (removeDOF)
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
    else if (removeDOF)
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
    std::cout << "Entropy in QH-approximation: " << entropy_sho << " cal / (mol * K)" << std::endl;
    std::cout << "Starting corrections for anharmonicity and Mutual Information to second order." << std::endl;

    //Corrections for anharmonicity and M.I.
    // I. Create PCA-Modes matrix
    Matrix_Class eigenvectors_t(transposed(eigenvectors));
    Matrix_Class pca_modes = eigenvectors_t * input;
    Matrix_Class entropy_anharmonic(pca_modes.rows(), 1u, 0.);
    Matrix_Class entropy_kNN(pca_modes.rows(), 1u, 0.);
    Matrix_Class entropy_mi(pca_modes.rows(), pca_modes.rows(), 0.);
    Matrix_Class statistical_entropy(pca_modes.rows(), 1u, 0.);
    Matrix_Class classical_entropy(pca_modes.rows(), 1u, 0.);


    // Modify PCA modes as the PCA eigenvalues have been modified. This is not detailed in the original paper
    // but sensible and reasonable to obtain valid values.
    pow(pca_modes, -1.);
    pca_modes = pca_modes * 1.05457172647 * 10e-34 / (sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp));
    //


    // II. Calculate Non-Paramteric Entropies
    // Marginal
#ifdef _OPENMP
#pragma omp parallel for firstprivate(pca_modes, k_kNN) shared(entropy_anharmonic, entropy_mi)
#endif
    for (int i = 0; i < (int)pca_modes.rows(); i++)
    {
      float_type* buffer = new float_type[k_kNN];

      KahanAccumulation<float_type> summedDistances;
      {
        for (size_t k = 0; k < pca_modes.cols(); k++)
        {
          float_type distance = log(sqrt(knn_distance_eucl_squared(pca_modes, 1, k_kNN, std::vector<size_t>{static_cast<size_t>(i)}, k, buffer)));
          summedDistances = KahanSum(summedDistances, distance);
        }
      }


      float_type distance = summedDistances.sum;

      distance /= float_type(pca_modes.cols()); // Number of draws
      float_type temp = float_type(pca_modes.cols()) * pow(3.14159265358979323846, 0.5);
      temp /= tgamma((1. / 2.) + 1.);
      distance += log(temp);
      temp = 0;
      if (k_kNN != 1u)
      {
        for (size_t l = 1; l < k_kNN; l++)
        {
          temp += 1.0 / float_type(l);
        }
      }
      distance -= temp;
      distance += 0.5772156649015328606065;
      entropy_kNN(i, 0u) = distance;
      distance = 0;

      //MI
      // Mutual Information Correction
      for (size_t j = i + 1; j < pca_modes.rows(); j++)
      {
        std::vector<size_t> query_rows{ (size_t)i,j };
        KahanAccumulation<float_type> kahan_acc;
        for (size_t k = 0; k < pca_modes.cols(); k++)
        {
          float_type dist_temp = log(sqrt(knn_distance_eucl_squared(pca_modes, 2, k_kNN, query_rows, k, buffer)));
          kahan_acc = KahanSum(kahan_acc, dist_temp);
        }
        distance = kahan_acc.sum;
        distance *= 2.;
        distance /= float_type(pca_modes.cols());
        float_type temp2 = float_type(pca_modes.cols()) * pow(3.14159265358979323846, 1);
        // Gamma(2) = 1. We omit temp2 /= 1.

        distance += log(temp2);
        temp2 = 0;
        if (k_kNN != 1)
        {
          for (size_t u = 1; u < k_kNN; u++)
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
    for (size_t i = 0; i < entropy_kNN.rows(); i++)
    {
      for (size_t j = (i + 1); j < entropy_kNN.rows(); j++)
      {
        if (
          pca_frequencies(i, 0u) <  ( temperatureInKelvin * 1.380648813 * 10e-23 / (1.05457172647 * 10e-34))
          && pca_frequencies(j, 0u) < ( temperatureInKelvin * 1.380648813 * 10e-23 / (1.05457172647 * 10e-34))
          )
        {
          entropy_mi(i, j) = entropy_kNN(i, 0u) + entropy_kNN(j, 0u) - entropy_mi(i, j);
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
          std::cout << "Notice: PCA-Modes " << i << " & " << j << " not corrected for M.I. since they are not in the classical limit\n";
          entropy_mi(i, j) = 0.0;
        }
      }
    }
    if (counterForLargeNegativeM_I_Terms > 0u)
    {
      std::cout << "Notice: Large negative M.I. term(s) detected. Check frequency of data sampling. (Do not worry, terms <0.0 are ignored anyway)\n";
    }
    for (size_t i = 0; i < entropy_kNN.rows(); i++)
    {
      statistical_entropy(i, 0u) = /*-1.0*  */ (log(alpha_i(i, 0u)) + log(sqrt(2. * 3.14159265358979323846 * 2.71828182845904523536)));
      classical_entropy(i, 0u) = -1.0 * (log(alpha_i(i, 0u)) - 1.);
      //statistical_entropy(i, 0u) = -1.0 * (log(pca_frequencies(i, 0u)) + log(sqrt(2 * 3.14159265358979323846 * 2.71828182845904523536)));
      //classical_entropy(i, 0u) = -1.0 * (log(pca_frequencies(i, 0u)) - 1.);

      entropy_anharmonic(i, 0u) = statistical_entropy(i, 0u) - entropy_kNN(i, 0u);

      // Debug output for developers
      if (Config::get().general.verbosity >= 5)
      {
        std::cout << "mode " << i << ". entropy kNN: " << entropy_kNN(i, 0u) << "\n";
        std::cout << "mode " << i << ". entropy anharmonic correction: " << entropy_anharmonic(i, 0u) << "\n";
        std::cout << "mode " << i << ". classical entropy: " << classical_entropy(i, 0u) << "\n";
        std::cout << "mode " << i << ". statistical entropy: " << statistical_entropy(i, 0u) << "\n";
        std::cout << "mode " << i << ". quantum entropy: " << quantum_entropy(i, 0u) << "\n";
        std::cout << "mode " << i << ". pca freq: " << pca_frequencies(i, 0u) << "\n";
        std::cout << "mode " << i << ". alpha (dimensionless, standard deviation): " << alpha_i(i, 0u) << "\n";
        std::cout << "mode " << i << ". standard deviation in mw-pca-units: " << sqrt(eigenvalues(i, 0u)) << "\n";
      }

      if (pca_frequencies(i, 0u) < (temperatureInKelvin * 1.380648813 * 10e-23 / (1.05457172647 * 10e-34)))
      {
        if (abs(entropy_anharmonic(i, 0u) / quantum_entropy(i, 0u)) < 0.007)
        {
          entropy_anharmonic(i, 0u) = 0.0;
          std::cout << "Notice: PCA-Mode " << i << " not corrected for anharmonicity (value too small)\n";
        }
      }
      else
      {
        std::cout << "Notice: PCA-Mode " << i << " not corrected for anharmonicity since it is not in the classical limit\n";
        entropy_anharmonic(i, 0u) = 0.0;
      }

      // Change dimensionless entropy to cal / K * mol
      entropy_anharmonic(i, 0u) *= 1.380648813 * 6.02214129 * 0.239005736;

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
}