
#include "TrajectoryMatrixClass.h"
//
//
using float_type = coords::float_type;
//
//
TrajectoryMatrixRepresentation::TrajectoryMatrixRepresentation(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords, \
  std::size_t start_frame_num, std::size_t offset, std::vector<std::size_t> trunc_atoms, std::size_t io_verbosity, \
  std::vector<std::size_t> internal_dih, int ref_frame_alignment, bool massweight)
{
  generateCoordinateMatrix(ci, coords,start_frame_num,offset,trunc_atoms,io_verbosity,internal_dih,ref_frame_alignment,massweight);
}

TrajectoryMatrixRepresentation::TrajectoryMatrixRepresentation(std::string const& filepath, std::size_t start_frame_num, std::size_t offset, std::vector<std::size_t> trunc_atoms)
{
  generateCoordinateMatrixfromPCAModesFile(filepath, start_frame_num,offset,trunc_atoms);
}

void TrajectoryMatrixRepresentation::generateCoordinateMatrixfromPCAModesFile(std::string const& filepath, \
  std::size_t start_frame_num, std::size_t offset, std::vector<std::size_t> trunc_atoms) //Config::get().entropy.entropy_trunc_atoms_num
{
  std::cout << "Reading snapshots/samples/trajectory from CAST-PCA-File \"" << filepath << "\".\n";
  std::cout << "Most I/O options are ignored as the data is taken from the file pretty much as-is, take care!" << std::endl;
  std::cout << "\"trunc_atoms_num\" option may be used to select PCA Modes that are kept.\n";
  std::cout << "No alignment is performed, we hope the PCA data was properly aligned.\n" << std::endl;
  std::cout << "Cartesian PCA modes are always assumed.\n" << std::endl;
  std::ifstream pcafile(filepath);

  if (pcafile.good())
  {
    pca::PrincipalComponentRepresentation pcadata = pca::PrincipalComponentRepresentation(filepath);
    Matrix_Class const& pcamodes = pcadata.getModes(); // "Trajectory in PCA - Modes following (columns are frames, rows are modes)"
    this->coordsMatrix = pcamodes;
    //
    //start_frame_num = Config::get().entropy.entropy_start_frame_num;
    //offset = Config::get().entropy.entropy_offset;
    //
    if (start_frame_num > 0u)
    {
      std::cout << "Starting with " << start_frame_num << "th frame from PCA-Modes.\n";
      this->coordsMatrix.shed_cols(0u, start_frame_num);
    }
    if (offset != 1u) // this still needs to be TESTED!!
    {
      std::size_t const numOriginalFrames = this->coordsMatrix.cols();
      std::size_t kept = 1u;
      for (std::size_t i = 1u; i < numOriginalFrames; ++i)
      {
        if (i % offset == 0u)
        {
          std::size_t temp1 = offset - 1u;
          this->coordsMatrix.shed_cols(kept, temp1);
          kept++;
        }
      }
    }
    if (trunc_atoms.size() > 0u) // this still needs to be tested.
    {
      std::vector<std::size_t> allDimensions;
      for (std::size_t i = 0u; i < this->coordsMatrix.rows(); ++i)
      {
        allDimensions.push_back(i);
      }
      for (std::size_t i = 0u; i < trunc_atoms.size(); ++i)
      {
        std::size_t const currentMode = trunc_atoms.at(i);
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

void TrajectoryMatrixRepresentation::generateCoordinateMatrix(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords,\
  std::size_t start_frame_num, std::size_t offset, std::vector<std::size_t> trunc_atoms, std::size_t io_verbosity, \
  std::vector<std::size_t> internal_dih, int ref_frame_alignment, bool massweight)
{
  const bool use_internal = internal_dih.size() != 0u;
  const bool trunc_cartesian = trunc_atoms.size() != 0u;
  if (use_internal && trunc_cartesian)
  {
    std::cout << "ERROR: You may not specify runcated cartesian modes and internal modes. CHeck your Inputfile! Aborting!\n";
    throw std::runtime_error("INPUT ERROR");
    return;
  }
  // Initialize the reference frame (for alignment etc)
  coords::Coordinates coords_ref(coords);
  auto holder = (*ci).PES()[ref_frame_alignment == -1 ? 0u : ref_frame_alignment].structure.cartesian;
  coords_ref.set_xyz(holder);



  //////////////
  //
  //  A L I G N M E N T
  //
  // if ref_Frame_alignment != -1: Alignment to this frame is performed
  if (ref_frame_alignment != -1 && use_internal)
  {
    std::cout << "Alignment is (in this case) redundant since internal coordinates are used. Alignment is skipped. Check your INPUTFILE please.\n";
    std::cout << "Continuing anyway...\n";
  }
  // Translational alignment of the reference frame
  if (ref_frame_alignment != -1)
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
  if (use_internal)
  {
    coordsMatrix = Matrix_Class((size_t) /* explicitly casting to round down */ ((ci->size() - start_frame_num) / offset), \
      internal_dih.size() * 2u);
  }
  // If truncated cartesians are desired, this section will
  // handle it.
  else if (trunc_cartesian)
  {
    coordsMatrix = Matrix_Class((size_t) /* explicitly casting to round down */ ((ci->size() - start_frame_num) / offset), \
      trunc_atoms.size() * 3u);
  }
  // This section is used if *all* cartesians of all atoms are used
  else
  {
    coordsMatrix = Matrix_Class((size_t) /* explicitly casting to round down */ ((ci->size() - Config::get().entropy.entropy_start_frame_num) / Config::get().entropy.entropy_offset), \
      coords.atoms().size() * 3u);
  }

  if (io_verbosity >= 3)
  {
    if (offset != 1u)
    {
      std::cout << "Only every " + std::to_string(offset) + "'th frame from the input trajectory is used." << std::endl;
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
    if (use_internal)
    {
      for (size_t i = start_frame_num; j < coordsMatrix.rows(); ++j, i += offset)
      {
        auto holder2 = ci->PES()[i].structure.cartesian;
        coords.set_xyz(holder2);
        coordsMatrix.set_row(j, ::matop::transformToOneline(coords, internal_dih, true));
      }
    }
    // This section if cartesian coordinates are used
    // They will *later* be massweighted
    else
    {
      for (size_t i = start_frame_num; j < coordsMatrix.rows(); ++j, i += offset)
      {
        auto holder2 = ci->PES()[i].structure.cartesian;
        coords.set_xyz(holder2);
        // Translational and rotational alignment
        if (ref_frame_alignment != -1)
        {
          // Alignes center of mass
          align::centerOfGeometryAlignment(coords);
          // Rotational alignment
          align::kabschAlignment(coords, coords_ref);
        }
        coordsMatrix.set_row(j, ::matop::transformToOneline(coords, trunc_atoms, false));
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
  if (!Config::get().entropy.entropy_use_internal && massweight)
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
  if (cov_determ = cov_matr.determ(), abs(cov_determ) < 1e-90)
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

  cov_matr *= (1.38064813 * /* 1e-23 J/K */ temperatureInKelvin * 2.718281828459 * 2.718281828459 / (1.054571726 /* * 10^-34 Js */ * 1.054571726 * 10e-45));
  cov_matr = cov_matr + Matrix_Class::identity(cov_matr.rows(), cov_matr.cols());
  float_type entropy_sho = cov_matr.determ();

  entropy_sho = log(entropy_sho) * 0.5 * 1.38064813 * 6.02214129 * 0.239;
  //This stems from: k_B * Avogadro * (to_calories) *0.5

  std::cout << "Entropy in QH-approximation: " << entropy_sho << " cal / (mol * K)" << std::endl;
  return entropy_sho;
}
