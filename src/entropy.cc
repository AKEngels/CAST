#include "entropy.h"

float_type ardakaniCorrection1D(float_type const& globMin, float_type const& globMax, float_type const& currentPoint, float_type const& NNdistance)
{
  if (currentPoint - NNdistance * 0.5 < globMin && currentPoint != globMin)
    return currentPoint + NNdistance * 0.5 - globMin;
  else if (currentPoint + NNdistance * 0.5 > globMax && currentPoint != globMax)
    return globMax - (currentPoint - NNdistance * 0.5);
  return NNdistance;
}

void scalePCACoordinatesForQuasiHarmonicTreatment(Matrix_Class& modes, float_type const& temperatureInK)
{
  pow(modes, -1.);
  modes = modes * 1.05457172647 * 10e-34 / (sqrt(1.380648813 * 10e-23 * temperatureInK));
  //
}

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
    float_type* distanceList = buffer;
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
    float_type* distanceList = buffer;
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

  TrajectoryMatrixRepresentation::TrajectoryMatrixRepresentation(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords)
  {
    generateCoordinateMatrix(ci, coords);
  }

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
      align::centerOfGeometryAlignment(coords_ref);
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

    if (Config::get().general.verbosity >= 3)
    {
      if (Config::get().entropy.entropy_offset != 1u)
      {
        std::cout << "Only every " + std::to_string(Config::get().entropy.entropy_offset) + "'th frame from the input trajectory is used." << std::endl;
      }
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
          coordsMatrix.set_row(j, ::matop::transformToOneline(coords, Config::get().entropy.entropy_internal_dih, true));
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
            align::centerOfGeometryAlignment(coords);
            // Rotational alignment
            align::kabschAlignment(coords, coords_ref);
          }
          coordsMatrix.set_row(j, ::matop::transformToOneline(coords, Config::get().entropy.entropy_trunc_atoms_num, false));
        }
      }
    }
    //std::cout << coordsMatrix << std::endl;

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

  float_type TrajectoryMatrixRepresentation::karplus() const
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

  float_type TrajectoryMatrixRepresentation::schlitter(float_type const temperatureInKelvin) const
  {
    std::cout << "\nCommencing entropy calculation:\nQuasi-Harmonic-Approx. according to Schlitter (see: doi:10.1016/0009-2614(93)89366-P)" << std::endl;
    Matrix_Class cov_matr = transposed(coordsMatrix);
    cov_matr = cov_matr - Matrix_Class(coordsMatrix.cols(), coordsMatrix.cols(), 1.0) * cov_matr / static_cast<float_type>(coordsMatrix.cols());
    cov_matr = transposed(cov_matr) * cov_matr;
    cov_matr = cov_matr / static_cast<float_type>(coordsMatrix.cols());

    cov_matr *= (1.38064813 * /* 10e-23 J/K */ temperatureInKelvin * 2.718281828459 * 2.718281828459 / (1.054571726 /* * 10^-34 Js */ * 1.054571726 * 10e-45));
    cov_matr = cov_matr + Matrix_Class::identity(cov_matr.rows(), cov_matr.cols());
    float_type entropy_sho = cov_matr.determ();

    entropy_sho = log(entropy_sho) * 0.5 * 1.38064813 * 6.02214129 * 0.239;
    //This stems from: k_B * Avogadro * (to_calories) *0.5

    std::cout << "Entropy in QH-approximation: " << entropy_sho << " cal / (mol * K)" << std::endl;
    return entropy_sho;
  }

}
