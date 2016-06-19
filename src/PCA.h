#pragma once
#include "matop.h"
namespace pca
{
	class PrincipalComponentRepresentation
	{
	public:
		PrincipalComponentRepresentation(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords);

		PrincipalComponentRepresentation(std::string const& filenameOfPCAModesFile);


		void writePCAModesFile(std::string const& filename = "pca_modes.dat");
		void writeHistogrammedProbabilityDensity(std::string const& filename = "pca_histogrammed.dat");
		void writeStocksDelta(std::string const& filename = "pca_stocksdelta.dat");
		void readEigenvectors(std::string const& filename = "pca_modes.dat");
		void readModes(std::string const& filename = "pca_modes2.dat");
		void generateCoordinateMatrix(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords);
		void generatePCAEigenvectorsFromCoordinates();
		void generatePCAModesFromPCAEigenvectorsAndCoordinates();
	private:
		Matrix_Class modes;
		Matrix_Class eigenvectors;
		Matrix_Class eigenvalues;
		Matrix_Class coordinatesMatrix;
	};
}