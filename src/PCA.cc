#include "PCA.h"
namespace pca
{
	void PrincipalComponentRepresentation::generateCoordinateMatrix(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords)
	{
		Matrix_Class matrix_aligned;

		//First, adjust number of truncated atoms to be used to zero, in case truncation is not to be used (duh)
		if (!Config::get().PCA.pca_trunc_atoms_bool) Config::set().PCA.pca_trunc_atoms_num = std::vector<size_t>();

		coords::Coordinates coords_ref(coords);
		auto holder = ci->PES()[Config::get().PCA.pca_ref_frame_num].structure.cartesian;
		coords_ref.set_xyz(holder);
		//Constructs two coordinate objects and sets reference frame according to INPUTFILE

		//Perform translational alignment for reference frame
		if (Config::get().PCA.pca_alignment && !Config::get().PCA.pca_use_internal)
		{
			::matop::align::centerOfMassAlignment(coords_ref);
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
					if (coords.atoms(Config::get().PCA.pca_internal_dih[i]).number() == 1u
						|| coords.atoms(coords.atoms(Config::get().PCA.pca_internal_dih[i]).ibond()).number() == 1u
						|| coords.atoms(coords.atoms(Config::get().PCA.pca_internal_dih[i]).iangle()).number() == 1u
						|| coords.atoms(coords.atoms(Config::get().PCA.pca_internal_dih[i]).idihedral()).number() == 1u)
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
				std::iota(std::begin(Config::set().PCA.pca_trunc_atoms_num), std::end(Config::set().PCA.pca_trunc_atoms_num), 1);
			}
			else if (Config::get().general.verbosity > 2U) std::cout << "Using truncated cartesian coordinates.\n";
			if (Config::get().PCA.pca_ignore_hydrogen)
			{
				if (Config::get().general.verbosity > 2U) std::cout << "Excluding hydrogen atoms.\n";
				int sizeOfVector = static_cast<int>(Config::get().PCA.pca_trunc_atoms_num.size());
				for (int i = 0u; i < sizeOfVector; i++)
				{
					if (coords.atoms().atom(Config::get().PCA.pca_trunc_atoms_num[i] - 1).number() == 1u)
					{
						Config::set().PCA.pca_trunc_atoms_num.erase(Config::set().PCA.pca_trunc_atoms_num.begin() + i);
						sizeOfVector--;
						i--;
					}
				}
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
						::matop::align::centerOfMassAlignment(coords); //Alignes center of mass
						::matop::align::kabschAlignment(coords, coords_ref); //Rotates
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
		if (Config::get().general.verbosity > 2U) std::cout << "Performing PCA transformation. This might take quite a while.\n";
		Matrix_Class cov_matr = (transposed(this->coordinatesMatrix));
		Matrix_Class ones(this->coordinatesMatrix.cols(), this->coordinatesMatrix.cols(), 1.0);
		cov_matr = cov_matr - ones * cov_matr / static_cast<float_type>(this->coordinatesMatrix.cols());
		cov_matr = transposed(cov_matr) * cov_matr;
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
		this->modes = transposed(eigenvectors) * this->coordinatesMatrix;
	}
	
	void PrincipalComponentRepresentation::readEigenvectors(std::string const& filename)
	{
		if (Config::get().general.verbosity > 2U) std::cout << "Reading PCA eigenvectors from file pca_modes.dat .\n";
		std::string iAmNotImportant_YouMayDiscardMe;
		//readEigenvectorsAndModes(eigenvectors, pca_modes, iAmNotImportant_YouMayDiscardMe);
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

	}

	void PrincipalComponentRepresentation::readModes(std::string const& filename)
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

		for (size_t i = 0u; i < eigenvectors.rows(); i++)
		{
			std::getline(pca_modes_stream, line);
		}
		std::getline(pca_modes_stream, line);
		std::getline(pca_modes_stream, line);
		std::getline(pca_modes_stream, line);
		std::getline(pca_modes_stream, line);
		std::getline(pca_modes_stream, line);
		this->modes.resize(stoi(line.substr(7, 10)), stoi(line.substr(20, 10)));

		for (size_t i = 0u; i < this->modes.rows(); i++)
		{
			std::getline(pca_modes_stream, line);
			size_t whitespace = 0u, lastWhitespace = 0u;
			for (size_t j = 0u; j < this->modes.cols(); j++)
			{
				lastWhitespace = whitespace;
				whitespace = line.find(" ", lastWhitespace + 1u);

				this->modes(i, j) = stod(line.substr(lastWhitespace, whitespace - lastWhitespace));
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

		if (Config::get().general.verbosity > 2U) std::cout << "Writing PCA modes, eigenvalues and eigenvectors to file.\n";

		double sum_of_all_variances = 0.0;
		for (size_t i = 0; i < eigenvalues.rows(); i++)
		{
			sum_of_all_variances += eigenvalues(i);
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
		std::vector<size_t> dimensionsToBeUsedInHistogramming = Config::get().PCA.pca_dimensions_for_histogramming;

		if (dimensionsToBeUsedInHistogramming.size() > 2)
		{
			std::cout << "More than 2 dimensions for histogramming is not yet supported!\n";
			std::cout << "Working with 2 dimensions";
			dimensionsToBeUsedInHistogramming.resize(2);
		}
		else
		{
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

			std::cout << "Starting: Histogramming...";
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
	}

	PrincipalComponentRepresentation::PrincipalComponentRepresentation(std::unique_ptr<coords::input::format>& ci, ::coords::Coordinates& coords)
	{
		generateCoordinateMatrix(ci, coords);
		this->generatePCAEigenvectorsFromCoordinates();
		this->generatePCAModesFromPCAEigenvectorsAndCoordinates();
	}

	PrincipalComponentRepresentation::PrincipalComponentRepresentation(std::string const& filenameOfPCAModesFile)
	{
		this->readModes(filenameOfPCAModesFile);
		this->readEigenvectors(filenameOfPCAModesFile);
	}

  void PrincipalComponentRepresentation::writeStocksDelta(std::string const& filename)
  {
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
}