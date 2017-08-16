#include "PCA.h"
namespace pca
{
	void PrincipalComponentRepresentation::generateCoordinateMatrix(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords)
	{
    // Tell 'em what we are gonna do:
    if (Config::get().general.verbosity > 2u)
      std::cout << "Generating new coordinate matrix from input trajectory according to options in the Inputfile." << std::endl;

		Matrix_Class matrix_aligned;

		//First, adjust number of truncated atoms to be used to zero, in case truncation is not to be used (duh)
		if (!Config::get().PCA.pca_trunc_atoms_bool) Config::set().PCA.pca_trunc_atoms_num = std::vector<size_t>();

		coords::Coordinates coords_ref(coords);
    if (Config::get().PCA.pca_ref_frame_num >= ci->size())
      throw std::runtime_error("Error in config-option pca_ref_frame_num: Number higher than total number of simulation frames.");
		auto holder = ci->PES()[Config::get().PCA.pca_ref_frame_num].structure.cartesian;
		coords_ref.set_xyz(holder);
		//Constructs two coordinate objects and sets reference frame according to INPUTFILE

		//Perform translational alignment for reference frame
		if (Config::get().PCA.pca_alignment && !Config::get().PCA.pca_use_internal)
		{
			align::centerOfMassAlignment(coords_ref);
		}

		// truncate internal coordinates "PCA.pca_internal_dih" 
		// or cartesians "PCA.pca_trunc_atoms" respectively
		// especially if hydrogens should be omitted
		if (Config::get().PCA.pca_use_internal)
		{
			if (Config::get().general.verbosity > 2U) std::cout << "Using dihedral angles.\n";
			if (Config::get().PCA.pca_ignore_hydrogen)
			{
				if (Config::get().general.verbosity > 2U) std::cout << "Excluding dihedrals involving hydrogen.\n";
				int sizeOfVector = static_cast<int>(Config::get().PCA.pca_internal_dih.size());
				for (int i = 0u; i < sizeOfVector; i++)
				{
          size_t ii = Config::get().PCA.pca_internal_dih[i];
          size_t ci = coords.atoms(ii).i_to_a();
          size_t ibondpartner = coords.atoms(ii).ibond();
          size_t ianglepartner = coords.atoms(ii).iangle();
          size_t idihpartner = coords.atoms(ii).idihedral();

          size_t cbondpartner = coords.atoms(ibondpartner).i_to_a();
          size_t canglepartner = coords.atoms(ianglepartner).i_to_a();
          size_t cdihpartner = coords.atoms(idihpartner).i_to_a();

          size_t num = coords.atoms(ci).number();
          size_t numbnd = coords.atoms(cbondpartner).number();
          size_t numang = coords.atoms(canglepartner).number();
          size_t numdih = coords.atoms(cdihpartner).number();

					if (num == 1u	|| numbnd == 1u || numang == 1u || numdih == 1u)
					{
						Config::set().PCA.pca_internal_dih.erase(Config::set().PCA.pca_internal_dih.begin() + i);
						sizeOfVector--;
						i--;
					}
				}
			}
			matrix_aligned = Matrix_Class((size_t) /* explicitly casting to round down */ (
				(ci->size() - Config::get().PCA.pca_start_frame_num) / Config::get().PCA.pca_offset),
				Config::get().PCA.pca_internal_dih.size() * 2u);
		}
		else
		{
			if (!Config::get().PCA.pca_trunc_atoms_bool)
			{
				if (Config::get().general.verbosity > 2U) std::cout << "Using all cartesian coordinates.\n";
				Config::set().PCA.pca_trunc_atoms_num = std::vector<size_t>(coords.atoms().size());
				// Fill with 0, 1, 2,..., .
				std::iota(std::begin(Config::set().PCA.pca_trunc_atoms_num), std::end(Config::set().PCA.pca_trunc_atoms_num), 0);
			}
			else if (Config::get().general.verbosity > 2U) std::cout << "Using truncated cartesian coordinates.\n";
			if (Config::get().PCA.pca_ignore_hydrogen)
			{
				if (Config::get().general.verbosity > 2U) std::cout << "Excluding hydrogen atoms.\n";
				int sizeOfVector = static_cast<int>(Config::get().PCA.pca_trunc_atoms_num.size());

				for (int i = 0u; i < sizeOfVector; i++)
				{
					if (coords.atoms().atom(Config::get().PCA.pca_trunc_atoms_num[i]).number() == 1u)
					{
						Config::set().PCA.pca_trunc_atoms_num.erase(Config::set().PCA.pca_trunc_atoms_num.begin() + i);
						sizeOfVector--;
						i--;
					}
				}

        if (coords.atoms().size() != Config::get().PCA.pca_trunc_atoms_num.size())
          Config::set().PCA.pca_trunc_atoms_bool = true;
			}
			matrix_aligned = Matrix_Class((size_t) /* explicitly casting to round down */ \
				((ci->size() - Config::get().PCA.pca_start_frame_num) / \
					Config::get().PCA.pca_offset), (Config::get().PCA.pca_trunc_atoms_num.size()) * 3u);
		}

		/* j counts the (truncated) matrix access, i the frames in ci */
		{
			size_t j = 0;
			if (Config::get().PCA.pca_use_internal)
			{
				for (size_t i = Config::get().PCA.pca_start_frame_num; j < matrix_aligned.rows(); ++j, i += Config::get().PCA.pca_offset)
				{
					auto holder2 = ci->PES()[i].structure.intern;
					coords.set_internal(holder2);
					matrix_aligned.row(j) = ::matop::transformToOneline(coords, Config::get().PCA.pca_internal_dih, true);
				}
			}
			else
			{
				for (size_t i = Config::get().PCA.pca_start_frame_num; j < matrix_aligned.rows(); ++j, i += Config::get().PCA.pca_offset)
				{
					auto holder2 = ci->PES()[i].structure.cartesian;
					coords.set_xyz(holder2);
					if (Config::get().PCA.pca_alignment)        //Translational and rotational alignment
					{
						align::centerOfMassAlignment(coords); //Alignes center of mass
						align::kabschAlignment(coords, coords_ref); //Rotates
					}
					matrix_aligned.row(j) = ::matop::transformToOneline(coords, Config::get().PCA.pca_trunc_atoms_num, false).row(0u);
				}
			}
		}

		transpose(matrix_aligned);

		if (!Config::get().PCA.pca_use_internal)
		{
			coords_ref.set_xyz(holder);
			if (!Config::get().PCA.pca_trunc_atoms_bool)
			{
				::matop::massweight(matrix_aligned, coords_ref, false);
			}
			else
			{
				::matop::massweight(matrix_aligned, coords_ref, false, Config::get().PCA.pca_trunc_atoms_num);
			}
		}
		this->coordinatesMatrix = matrix_aligned;
	}

	void PrincipalComponentRepresentation::generatePCAEigenvectorsFromCoordinates()
	{
    // Tell 'em what we are gonna do:
    if (Config::get().general.verbosity > 2u)
      std::cout << "Generating PCA Eigenvectors from coordinate matrix." << std::endl;

		if (Config::get().general.verbosity > 2U) std::cout << "Performing PCA transformation. This might take quite a while.\n";
		Matrix_Class cov_matr = (transposed(this->coordinatesMatrix));
		Matrix_Class ones(this->coordinatesMatrix.cols(), this->coordinatesMatrix.cols(), 1.0);
		cov_matr = Matrix_Class (cov_matr - ones * cov_matr / static_cast<float_type>(this->coordinatesMatrix.cols()));
		cov_matr = Matrix_Class (transposed(cov_matr) * cov_matr);
		cov_matr = cov_matr / static_cast<float_type>(this->coordinatesMatrix.cols());
		float_type cov_determ = 0.;
		int *cov_rank = new int;
		cov_matr.eigensym(this->eigenvalues, this->eigenvectors, cov_rank);
		if (*cov_rank < (int)eigenvalues.rows() || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90))
		{
			std::cout << "Notice: covariance matrix is singular.\n";
			std::cout << "Details: rank of covariance matrix is " << *cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << ".\n";
		}
		delete cov_rank;
	}

	void PrincipalComponentRepresentation::generatePCAModesFromPCAEigenvectorsAndCoordinates()
	{
    // Tell 'em what we are gonna do:
    if (Config::get().general.verbosity > 2u)
      std::cout << "Generating PCA modes from coordinate matrix and PCA Eigenvectors." << std::endl;

		this->modes = Matrix_Class (transposed(eigenvectors) * this->coordinatesMatrix);
	}
	
	void PrincipalComponentRepresentation::readEigenvectors(std::string const& filename)
	{
		if (Config::get().general.verbosity > 2U) std::cout << "Reading PCA eigenvectors from file" << filename << "." << std::endl;
		std::string iAmNotImportant_YouMayDiscardMe;
		std::ifstream pca_modes_stream(filename, std::ios::in);
		std::string line;
		std::getline(pca_modes_stream, line);
		while (line.find("Eigenvectors") == std::string::npos)
		{
			std::getline(pca_modes_stream, line);
		}

		std::getline(pca_modes_stream, line);

    eigenvectors.resize(stoi(line.substr(7, 16)), stoi(line.substr(20, 10)));

    for (size_t i = 0u; i < eigenvectors.rows(); i++)
    {
      std::getline(pca_modes_stream, line);
      for (size_t j = 0u; j < eigenvectors.cols(); j++)
      {
        std::string number = line.substr(j * 19u, 19u);
        eigenvectors(i, j) = stod(number);
      }
    }

	}

	void PrincipalComponentRepresentation::readModes(std::string const& filename)
	{
    if (Config::get().general.verbosity > 2U) std::cout << "Reading PCA modes from file" << filename << "." << std::endl;
		std::ifstream pca_modes_stream(filename, std::ios::in);
		std::string line;
		std::getline(pca_modes_stream, line);
		while (line.find("Eigenvectors") == std::string::npos)
		{
			std::getline(pca_modes_stream, line);
		}

		std::getline(pca_modes_stream, line);

		for (size_t i = 0u; i < eigenvectors.rows(); i++)
		{
			std::getline(pca_modes_stream, line);
		}
		std::getline(pca_modes_stream, line);
		std::getline(pca_modes_stream, line);
		std::getline(pca_modes_stream, line);
		std::getline(pca_modes_stream, line);
		std::getline(pca_modes_stream, line);
    modes.resize(stoi(line.substr(7, 16)), stoi(line.substr(20, 10)));

    for (size_t i = 0u; i < modes.rows(); i++)
    {
      std::getline(pca_modes_stream, line);
      for (size_t j = 0u; j < modes.cols(); j++)
      {
        std::string number = line.substr(j * 19u, 19u);
        modes(i, j) = stod(number);
      }
    }
	}

	void PrincipalComponentRepresentation::writePCAModesFile(std::string const& filename)
	{
		///////////////////////////////////////
		// Build the additional information string which is necessary for string i/o in PCAproc task
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

		if (Config::get().general.verbosity > 2U) std::cout << "Writing PCA modes, eigenvalues and eigenvectors to file " << filename << "." << std::endl;

		double sum_of_all_variances = 0.0;
		for (size_t i = 0; i < eigenvalues.rows(); i++)
		{
			sum_of_all_variances += eigenvalues(i, 0u);
		}

		std::cout << "Working with " << modes.rows() << " dimensions.\n";

		std::ofstream pca_modes_stream(filename, std::ios::out);
		pca_modes_stream << "Working with " << modes.rows() << " dimensions.\n";
		pca_modes_stream << "PCA-Modes are ordered with descending variances.\n\n\n\n";
		pca_modes_stream << "Eigenvalues of Covariance Matrix:\n";
		pca_modes_stream << eigenvalues << "\n\n\n";
		pca_modes_stream << "Eigenvectors of Covariance Matrix:\nSize = " << std::setw(10) << eigenvectors.rows() << " x " << std::setw(10) << eigenvectors.cols() << "\n";
		pca_modes_stream << eigenvectors << "\n\n\n";
		pca_modes_stream << "Trajectory in PCA - Modes following (columns are frames, rows are modes)\nSize = " << std::setw(10) << modes.rows() << " x " << std::setw(10) << modes.cols() << "\n";
		pca_modes_stream << modes << "\n";
		pca_modes_stream << "Additional Information following:\n" << additionalInformation << "\n";
	}

	void PrincipalComponentRepresentation::writeHistogrammedProbabilityDensity(std::string const& filename)
	{
		using namespace histo;
    if (!Config::get().PCA.pca_histogram_all_marginal_degrees_of_freedom)
    {
      std::vector<size_t> dimensionsToBeUsedInHistogramming = Config::get().PCA.pca_dimensions_for_histogramming;

      DimensionalHistogram<float_type> *histograms_p;
      if (Config::get().PCA.pca_histogram_number_of_bins > 0u)
      {
        size_t histogramBins = Config::get().PCA.pca_histogram_number_of_bins;
        //histograms_p = new DimensionalHistogram<float_type>((size_t)pca_modes.rows(), histogramBins);
        histograms_p = new DimensionalHistogram<float_type>((size_t)dimensionsToBeUsedInHistogramming.size(), histogramBins);
      }
      else if (Config::get().PCA.pca_histogram_width > 0.)
      {
        float_type histogramWidth = Config::get().PCA.pca_histogram_width;
        histograms_p = new DimensionalHistogram<float_type>((size_t)dimensionsToBeUsedInHistogramming.size(), histogramWidth);
      }
      else
      {
        throw "Error in output of probability density, exiting.\n You need to specify either a number of histogram bins or a bin width in the inputfile.\n";
      }

      if (Config::get().general.verbosity > 2U)
        std::cout << "Writing histogrammed probability density to file " << filename << "." << std::endl;

      //Filling the histogram
      for (size_t j = 0u; j < modes.cols(); j++)
      {
        std::vector<float_type> temp(dimensionsToBeUsedInHistogramming.size());
        for (size_t i = 0u; i < dimensionsToBeUsedInHistogramming.size(); i++)
        {
          temp[i] = modes(dimensionsToBeUsedInHistogramming[i] - 1, j);
        }
        histograms_p->add_value(temp);
      }

      histograms_p->distribute();
      histograms_p->writeProbabilityDensity(filename);
      histograms_p->writeAuxilaryData("auxdata_" + filename);
      delete histograms_p;
    }
    else
    {
      std::vector<size_t> dimensionsToBeUsedInHistogramming = { 1u };

      if (Config::get().general.verbosity > 2U)
        std::cout << "Writing histogrammed probability density to file " << filename << "." << std::endl;
      //Filling the histograms
      for (size_t currentDimension = 0u; currentDimension < modes.rows(); currentDimension++)
      {
        DimensionalHistogram<float_type> *histograms_p;
        if (Config::get().PCA.pca_histogram_number_of_bins > 0u)
        {
          size_t histogramBins = Config::get().PCA.pca_histogram_number_of_bins;
          //histograms_p = new DimensionalHistogram<float_type>((size_t)pca_modes.rows(), histogramBins);
          histograms_p = new DimensionalHistogram<float_type>(1u, histogramBins);
        }
        else if (Config::get().PCA.pca_histogram_width > 0.)
        {
          float_type histogramWidth = Config::get().PCA.pca_histogram_width;
          histograms_p = new DimensionalHistogram<float_type>(1u, histogramWidth);
        }
        else
        {
          throw "Error in output of probability density, exiting.\n You need to specify either a number of histogram bins or a bin width in the inputfile.\n";
        }

        for (size_t j = 0u; j < modes.cols(); j++)
        {
          histograms_p->add_value(std::vector<float_type>{modes(currentDimension, j)});
        }
        histograms_p->distribute();
        histograms_p->writeProbabilityDensity("dim" + std::to_string(currentDimension + 1u) + "_" + filename);
        histograms_p->writeAuxilaryData("auxdata_dim" + std::to_string(currentDimension + 1u) + "_" + filename);
        delete histograms_p;
      }
    }
	}

	PrincipalComponentRepresentation::PrincipalComponentRepresentation(std::unique_ptr<coords::input::format>& ci, ::coords::Coordinates& coords)
	{
		generateCoordinateMatrix(ci, coords);
		this->generatePCAEigenvectorsFromCoordinates();
		this->generatePCAModesFromPCAEigenvectorsAndCoordinates();
	}

	PrincipalComponentRepresentation::PrincipalComponentRepresentation(std::string const& filenameOfPCAModesFile)
	{
    this->readEigenvectors(filenameOfPCAModesFile);
		this->readModes(filenameOfPCAModesFile);
	}

  void PrincipalComponentRepresentation::writeStocksDelta(std::string const& filename)
  {
    if (Config::get().general.verbosity > 2U) 
      std::cout << "Writing Stock's delta to file " << filename << "." << std::endl;


    std::ofstream stream(filename, std::ios::out);
    stream << "Stock's Delta; see DOI 10.1063/1.2746330\n\n";
    if (Config::get().PCA.pca_use_internal)
    {
      for (size_t i = 0u; i < this->modes.rows(); i++)
      {
        stream << "PCA Mode " << i << ":\n";
        for (size_t j = 0u; j < 2u * Config::get().PCA.pca_internal_dih.size(); j += 2u)
        {
          stream << "Dihedral " << std::setw(6) << Config::get().PCA.pca_internal_dih[j / 2u] << std::setw(0) << ": ";

          float_type delta =
            this->eigenvectors(i, j) * this->eigenvectors(i, j) +
            this->eigenvectors(i, j + 1) * this->eigenvectors(i, j + 1);

          stream << std::setw(10) << delta << std::setw(0) << "\n";
        }
        stream << "------------------------------\n";
      }
    }
    else
    {
      for (size_t i = 0u; i < this->modes.rows(); i++)
      {
        stream << "PCA Mode " << i << ":\n";
        for (size_t j = 0u; j < 3u * Config::get().PCA.pca_trunc_atoms_num.size(); j+=3u)
        {
          stream << "Atom " << std::setw(6) << Config::get().PCA.pca_trunc_atoms_num[j/3u] << std::setw(0) << ": ";

          float_type delta = 
            this->eigenvectors(i, j) * this->eigenvectors(i, j) + 
            this->eigenvectors(i, j + 1) * this->eigenvectors(i, j + 1) + 
            this->eigenvectors(i, j + 2) * this->eigenvectors(i, j + 2);

          stream << std::setw(10) << delta << std::setw(0) << "\n";
        }
        stream << "------------------------------\n";
      }
    }
  }

  ProcessedPrincipalComponentRepresentation::ProcessedPrincipalComponentRepresentation(std::string const& filenameOfPCAModesFile)
    : PrincipalComponentRepresentation(filenameOfPCAModesFile)
  {
    this->readAdditionalInformation(filenameOfPCAModesFile);
    this->restoreCoordinatesMatrix();
  }

  void ProcessedPrincipalComponentRepresentation::readAdditionalInformation(std::string const& filename)
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
    size_t vec_rows = static_cast<size_t>(stoi(line.substr(7, 16)));
    size_t vec_cols = static_cast<size_t>(stoi(line.substr(20, 10)));

    //eigenvectors.resize(stoi(line.substr(7, 16)), stoi(line.substr(20, 10)));

    for (size_t i = 0u; i < vec_rows; i++)
    {
      std::getline(pca_modes_stream, line);
    }
    std::getline(pca_modes_stream, line);
    std::getline(pca_modes_stream, line);
    std::getline(pca_modes_stream, line);
    std::getline(pca_modes_stream, line);
    std::getline(pca_modes_stream, line);

    //trajectory.resize(stoi(line.substr(7, 16)), stoi(line.substr(20, 10)));
    size_t traj_rows = static_cast<size_t>(stoi(line.substr(7, 16)));
    size_t traj_cols = static_cast<size_t>(stoi(line.substr(20, 10)));


    for (size_t i = 0u; i < traj_rows; i++)
    {
      std::getline(pca_modes_stream, line);
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
        std::cout << "Could not read additional Information from pca_modes file.\n";
      }
    }
  }

  void ProcessedPrincipalComponentRepresentation::determineStructures(std::unique_ptr<coords::input::format>& ci, ::coords::Coordinates& coords)
  {
    if (Config::get().PCA.proc_desired_start.size() > modes.rows() || Config::get().PCA.proc_desired_stop.size() > modes.rows())
    {
      std::cout << "Desired PCA-Ranges have higher dimensionality then modes. Omitting the last values.\n";
    }

    for (size_t j = 0u; j < coordinatesMatrix.cols(); j++)
    {
      bool isStructureInRange = true;
      for (size_t i = 0u; i < modes.rows() && i < std::max(Config::get().PCA.proc_desired_stop.size(), Config::get().PCA.proc_desired_start.size()); i++)

      {
        if (i < Config::get().PCA.proc_desired_start.size())
        {
          if (modes(i, j) < Config::get().PCA.proc_desired_start[i])
          {
            isStructureInRange = false;
          }
        }
        if (i < Config::get().PCA.proc_desired_stop.size())
        {
          if (modes(i, j) > Config::get().PCA.proc_desired_stop[i])
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

    if (Config::get().general.verbosity > 1u) std::cout << "Found " << structuresToBeWrittenToFile.size() << " structures in desired range.\nWorking..." << std::endl;

    // Case: Cartesian Coordinates.
    if (additionalInformation.substr(0, 3) == "car")
    {
      //Additional Information Processing -> Read "DOFS that were used" from file. Put their identifying numbers in vector.
      std::stringstream ss(additionalInformation.substr(4));
      std::string buffer;
      std::vector<size_t> tokens;
      std::deque<bool> alreadyFoundStructures(ci->size(), false);

      while (ss >> buffer) tokens.push_back((size_t)std::stoi(buffer));

      if (tokens.size() == 0u)
      {
        tokens.resize(coords.atoms().size());
        std::iota(std::begin(tokens), std::end(tokens), 0);
      }

      if (tokens.size() != 0u)
      {
        ::matop::undoMassweight(coordinatesMatrix, coords, false, tokens);
        for (size_t i = 0u; i < structuresToBeWrittenToFile.size(); i++)
        {
          Matrix_Class out_mat(3, coordinatesMatrix.rows() / 3u);
          for (size_t j = 0u; j < coordinatesMatrix.rows(); j = j + 3)
          {
            out_mat(0, j / 3u) = coordinatesMatrix(j, structuresToBeWrittenToFile[i]);
            out_mat(1, j / 3u) = coordinatesMatrix(j + 1u, structuresToBeWrittenToFile[i]);
            out_mat(2, j / 3u) = coordinatesMatrix(j + 2u, structuresToBeWrittenToFile[i]);
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
            this->foundStructures.push_back(ci->PES()[structureNumber]);
            alreadyFoundStructures[structureNumber] = true;
          }
          else
          {
            std::cout << "Could not find structure restored from PCA-Modes in ensemble of structures from original coordinates.\n";
            std::cout << "This means that there was no provided structure with a deviance of less than 0.5% to the current restored structure.\n\n";
          }
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
              float_type compareFromPCA1 = coordinatesMatrix(quicksearch, structuresToBeWrittenToFile[i]);
              float_type compareFromTrajectory1 = std::cos(coords.intern(j).azimuth().radians());
              float_type compareFromPCA2 = coordinatesMatrix(quicksearch + 1u, structuresToBeWrittenToFile[i]);
              float_type compareFromTrajectory2 = std::sin(coords.intern(j).azimuth().radians());
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
          //If match was found, store it.
          if (structureFound)
          {
            alreadyFoundStructures[structureNumber] = true;
            //coords.set_pes(ci->PES()[k]);
            foundStructures.push_back(ci->PES()[k]);
            //outstream << coords;
          }
        }
        if (!structureFound)
        {
          throw std::runtime_error(
            "Could not find structure restored from PCA-Modes in ensemble of structures from original coordinates.\nYou probably made a mistake somewhere in your INPUTFILE.\nDo not consider structures written out after this message as valid.\n"
            );
        }
      }
    }
  }

  void ProcessedPrincipalComponentRepresentation::restoreCoordinatesMatrix()
  {
    //Undoing PCA
    this->coordinatesMatrix = Matrix_Class (eigenvectors * modes);
  }

#ifdef CAST_USE_ARMADILLO
  void PrincipalComponentRepresentation::createGaussianMixtureModel(std::vector<size_t> const& dimensions, size_t numberOfGaussians)
  {
    using namespace arma;
    Matrix_Class submodes(dimensions.size(), this->modes.cols());
    for (unsigned int i = 0u; i < dimensions.size(); i++)
    {
      for (unsigned int j = 0u; j < this->modes.cols(); j++)
      {
        submodes(i, j) = this->modes(dimensions.at(i), j);
      }
    }

    arma::Mat<double> & pca_modes_matrix = submodes;
    arma::gmm_full model;

    bool status = model.learn(pca_modes_matrix, numberOfGaussians, eucl_dist, random_spread, 10, 7, 0, true);
    if (status == false)
    {
      std::cout << "learning failed" << std::endl;
    }

    Matrix_Class means = model.means;
    for (unsigned int i = 0u; i < means.cols(); i++)
    {
      std::cout << "Mean of Gaussian " << i << ":\n";
      for (unsigned int j = 0u; j < means.rows(); j++)
      {
        std::cout << means(j, i) << " ";
      }
      std::cout << "\n";
    }

    for (unsigned int i = 0u; i < means.cols(); i++)
    {
      std::cout << "Covariance Matrix of Gaussian " << i << ":\n";
      Matrix_Class current = model.fcovs.slice(i);
      std::cout << current << "\n";
    }
    std::cout << "\n";
    for (unsigned int i = 0u; i < means.cols(); i++)
    {
      std::cout << "Weight of Gaussian " << i << ": ";
      std::cout << model.hefts(i) << "\n";
    }
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";
  }
#endif

  void ProcessedPrincipalComponentRepresentation::writeDeterminedStructures(::coords::Coordinates const& coord_in, std::string const& filenameExtension)
  {
    std::ofstream outstream(coords::output::filename(filenameExtension).c_str(), std::ios::app);
    coords::Coordinates coords(coord_in);

    for (size_t i = 0u; i < this->foundStructures.size(); ++i)
    {
      coords.set_pes(this->foundStructures[i]);
      coords.to_internal();
      outstream << coords;
    }
  }
}