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

  coords::float_type calculateRotiationalEntropy(coords::Coordinates const& coordobj, coords::float_type temperatureInK, std::size_t symmetryNumber)
  {
    using coords::float_type;
    float_type rotEntropy = 1.;
    Matrix_Class inertiaTensor = coordobj.getInertiaTensor();
    std::cout << "\n" << inertiaTensor << "\n";
    std::pair< Matrix_Class, Matrix_Class> eigReturn = inertiaTensor.eigensym(true,true);
    Matrix_Class const& eigenval = eigReturn.first;
    std::cout << "\n" << eigReturn.first << "\n";
    std::cout << "\n" << eigReturn.second << "\n";
    eigReturn.first = Matrix_Class(3u,1u);
    eigReturn.first(0,0) = 307.416;
    eigReturn.first(1,0) = 50.54344;
    //eigReturn.first(1,0) = 1050.54344;
    eigReturn.first(2, 0) = 35.36481;
    //eigReturn.first(2,0) = 1335.36481;
    temperatureInK = 298.15;
    // via http://www.codessa-pro.com/descriptors/thermodynamic/rem.htm
    double inertiaProduct = 1.;
    for (std::size_t i = 0u; i < 3u; ++i)
    {
      std::cout << eigenval(i,0u) << " " << std::endl;
      const double inertia = eigenval(i,0u) / 1000.0 /* to kg/(mol*A²) */ / constants::N_avogadro /* to kg/A² */ / (constants::angstrom2meters * constants::angstrom2meters);
      inertiaProduct *= inertia;
      
    }
    rotEntropy = std::sqrt(inertiaProduct);
    rotEntropy = std::sqrt(constants::pi) / static_cast<coords::float_type>(symmetryNumber) * std::pow(temperatureInK,3.0/2.0) / rotEntropy;
    rotEntropy = std::log(rotEntropy) + 3.0/2.0;
    rotEntropy *= constants::gas_constant_R_kcal_per_mol_kelvin;

    return rotEntropy;
  }

  coords::float_type calculateTranslationalEntropy(coords::float_type totalMass, coords::float_type temperatureInK, coords::float_type volumeOfSystem)
  {
    using coords::float_type;
    constexpr float_type referencePressure = 101325.0; // Pascal
    float_type transEntropy = 1.;
    totalMass /= constants::N_avogadro;
    totalMass /= 1000.0; // mass in kg
    transEntropy *= std::pow(2*constants::pi* totalMass* constants::boltzmann_constant_kb_SI_units * temperatureInK /(constants::h_SI_units* constants::h_SI_units),1.5);
    transEntropy *= constants::boltzmann_constant_kb_SI_units * temperatureInK / referencePressure;
    transEntropy = std::log(transEntropy);
    transEntropy += 5.0 / 2.0;
    transEntropy *= constants::gas_constant_R_kcal_per_mol_kelvin;
    return transEntropy;
  }

  coords::float_type calculateTranslationalEntropy(coords::Coordinates const& coords, coords::float_type temperatureInK, coords::float_type volumeOfSystem)
  {
    float_type mass = 0.;
    for (auto const& i : coords.atoms())
    {
      mass += i.mass();
    }
    return calculateTranslationalEntropy(mass,temperatureInK,volumeOfSystem);
  }

  TrajectoryMatrixRepresentation::TrajectoryMatrixRepresentation(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords)
  {
    generateCoordinateMatrix(ci, coords);
  }

  TrajectoryMatrixRepresentation::TrajectoryMatrixRepresentation(std::string const& filepath)
  {
    generateCoordinateMatrixfromPCAModesFile(filepath);
  }

  void TrajectoryMatrixRepresentation::generateCoordinateMatrixfromPCAModesFile(std::string const& filepath)
  {
    std::cout << "Reading snapshots/samples/trajectory from CAST-PCA-File \"" << filepath << "\".\n";
    std::cout << "Most I/O options are ignored as the data is taken from the file pretty much as-is, take care!" << std::endl;
    std::cout << "\"trunc_atoms_num\" option is used to select PCA Modes that are kept.\n";
    std::cout << "No alignment is performed, we hope the PCA data was properly aligned.\n" << std::endl;
    std::cout << "Cartesian PCA modes are always assumed.\n" << std::endl;
    std::ifstream pcafile(filepath);
    
    if (pcafile.good())
    {
      pca::PrincipalComponentRepresentation pcadata = pca::PrincipalComponentRepresentation(filepath); 
      Matrix_Class const& pcamodes = pcadata.getModes(); // "Trajectory in PCA - Modes following (columns are frames, rows are modes)"
      this->coordsMatrix = pcamodes;
      //Config::get().entropy.entropy_start_frame_num
      if(Config::get().entropy.entropy_start_frame_num > 0u)
      {
        std::cout << "Starting with " << Config::get().entropy.entropy_start_frame_num << "th frame from PCA-Modes.\n";
        this->coordsMatrix.shed_cols(0u, Config::get().entropy.entropy_start_frame_num);
      }
      if (Config::get().entropy.entropy_offset != 1u) // this still needs to be TESTED!!
      {
        std::size_t const numOriginalFrames = this->coordsMatrix.cols();
        std::size_t kept = 1u;
        for (std::size_t i = 1u; i < numOriginalFrames; ++i)
        {
          if (i % Config::get().entropy.entropy_offset == 0u)
          {
            std::size_t temp1 = Config::get().entropy.entropy_offset - 1u;
            this->coordsMatrix.shed_cols(kept, temp1);
            kept++;
          }
        }
      }
      if (Config::get().entropy.entropy_trunc_atoms_num.size() != 0u && Config::get().entropy.entropy_trunc_atoms_bool) // this still needs to be tested.
      {
        std::vector<std::size_t> allDimensions;
        for (std::size_t i = 0u; i < this->coordsMatrix.rows(); ++i)
        {
          allDimensions.push_back(i);
        }
        for (std::size_t i = 0u; i < Config::get().entropy.entropy_trunc_atoms_num.size(); ++i)
        {
          std::size_t const currentMode = Config::get().entropy.entropy_trunc_atoms_num.at(i);
          if (std::find(allDimensions.begin(), allDimensions.end(), currentMode) != allDimensions.end()) 
          {
            // via https://stackoverflow.com/questions/3385229/c-erase-vector-element-by-value-rather-than-by-position
            allDimensions.erase(std::remove(allDimensions.begin(), allDimensions.end(), currentMode), allDimensions.end());
          }
          else 
          {
            throw std::runtime_error("Entry specified in input-file option \"entropy_trunc_atoms_num\" is out of bounds.");
          }
          // Sort in descending order
        }
        // via https://stackoverflow.com/questions/9025084/sorting-a-vector-in-descending-order
        std::sort(allDimensions.rbegin(), allDimensions.rend());
        for (std::size_t i = 0u; i < allDimensions.size(); ++i)
        {
          this->coordsMatrix.shed_row(allDimensions.at(i));
        }
      }
    }
    else
    {
      throw std::runtime_error("Path to pca-modes-file not valid. Aborting.");
    }

    // Rows are DOFs, Columns are frames!

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
