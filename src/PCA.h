#pragma once
#include "Matrix_Class.h"
#include "histogram.h"
#include "alignment.h"
namespace pca
{
  using Matrix_Class = ::Matrix_Class;
  class PrincipalComponentRepresentation
  {
  public:

    //using float_type = coords::float_type;

    /**default constructor
    after using this you have to create eigenvectors and modes, either by calculating them or by reading from file*/
    PrincipalComponentRepresentation() {};
    /**constructor that calculates eigenvectors and modes from coordinates*/
    PrincipalComponentRepresentation(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords);
    /**constructor that reads eigenvectors and modes from file (normally 'pca_modes.dat')
    you need to create the coordinate matrix seperately for printing stock's delta*/
    PrincipalComponentRepresentation(std::string const& filenameOfPCAModesFile);

    void writePCAModesFile(std::string const& filename = "pca_modes.dat");
    void writeHistogrammedProbabilityDensity(std::string const& filename = "pca_histogrammed.dat");
    void writeStocksDelta(std::string const& filename = "pca_stocksdelta.dat");
    void readEigenvectors(std::string const& filename = "pca_modes.dat");
    void readModes(std::string const& filename = "pca_modes.dat");
    void readEigenvalues(std::string const& filename = "pca_modes.dat");
    void generateCoordinateMatrix(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords);
    void generatePCAEigenvectorsFromCoordinates();
    void generatePCAModesFromPCAEigenvectorsAndCoordinates();
    Matrix_Class const& getModes() const;
  protected:
    Matrix_Class modes;
    Matrix_Class eigenvectors;
    Matrix_Class eigenvalues;
    Matrix_Class mw_coordinatesMatrix;
  };

  class ProcessedPrincipalComponentRepresentation
    : public PrincipalComponentRepresentation
  {
  public:
    ProcessedPrincipalComponentRepresentation(std::string const& filenameOfPCAModesFile);
    void readAdditionalInformation(std::string const& filename = "pca_modes.dat");
    void determineStructures(std::unique_ptr<coords::input::format>& ci, ::coords::Coordinates& coords);
    void restoreCoordinatesMatrix();
    void writeDeterminedStructures(::coords::Coordinates const& coord_in, std::string const& filenameExtension = "_pca_selection");
  private:
    std::vector<size_t> structuresToBeWrittenToFile;
    std::string additionalInformation;
    std::vector<coords::PES_Point> foundStructures;
  };
}