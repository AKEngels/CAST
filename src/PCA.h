#pragma once
#include "Matrix_Class.h"
#include "histogram.h"
#include "alignment.h"
#include "scon/scon_serialization.h"
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

    void writePCAModesFile(std::string const& filename = "pca_modes.dat") const;
    void writePCAModesBinaryFile(std::string const& filename = "pca_modes.cbf") const;
    void writeHistogrammedProbabilityDensity(std::string const& filename = "pca_histogrammed.dat");
    void writeStocksDelta(std::string const& filename = "pca_stocksdelta.dat");
    void readEigenvectors(std::string const& filename = "pca_modes.dat");
    void readModes(std::string const& filename = "pca_modes.dat");
    void readEigenvalues(std::string const& filename = "pca_modes.dat");
    void readBinary(std::string const& filename, bool readEigenvalues = true, bool readEigenvectors = true, bool readTrajectory = true, bool readPCAModes = true);
    void generateCoordinateMatrix(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords);
    void generatePCAEigenvectorsFromCoordinates();
    void generatePCAModesFromPCAEigenvectorsAndCoordinates();
    Matrix_Class const& getModes() const;
    Matrix_Class const& getEigenvectors() const;
    Matrix_Class const& getEigenvalues() const;
    Matrix_Class const& getMWTrajectoryMatrix() const;
  protected:
    Matrix_Class modes;
    Matrix_Class eigenvectors;
    Matrix_Class eigenvalues;
    Matrix_Class mw_coordinatesMatrix;
  public:
    template<class Strm>
    friend scon::binary_stream<Strm>& operator<< (scon::binary_stream<Strm>& strm, PrincipalComponentRepresentation const& pca)
    {
      const std::vector<std::vector<float_type>> vec_modes = pca.getModes().to_std_vector();
      std::vector<std::vector<float_type>> const vec_eigenvectors = pca.getEigenvectors().to_std_vector();
      std::vector<std::vector<float_type>> const vec_eigenvalues = pca.getEigenvalues().to_std_vector();
      std::vector<std::vector<float_type>> const vec_mw_coordinatesMatrix = pca.getMWTrajectoryMatrix().to_std_vector();
      //
      std::array<std::size_t, 8u> const sizes = {
        vec_modes.size(), vec_modes.at(0u).size(), 
        vec_eigenvectors.size(), vec_eigenvectors.at(0u).size(),
        vec_eigenvalues.size(), vec_eigenvalues.at(0u).size(),
        vec_mw_coordinatesMatrix.size(), vec_mw_coordinatesMatrix.at(0u).size(),
        };
      //
      for (auto const& s : sizes)
        strm << s;
      for (std::size_t i = 0u; i < vec_modes.size(); ++i)
      {
        for (std::size_t j = 0u; j < vec_modes.at(0u).size(); ++j)
        {
          strm << vec_modes.at(i).at(j);
        }
      }
      for (std::size_t i = 0u; i < vec_eigenvectors.size(); ++i)
      {
        for (std::size_t j = 0u; j < vec_eigenvectors.at(0u).size(); ++j)
        {
          strm << vec_eigenvectors.at(i).at(j);
        }
      }
      for (std::size_t i = 0u; i < vec_eigenvalues.size(); ++i)
      {
        for (std::size_t j = 0u; j < vec_eigenvalues.at(0u).size(); ++j)
        {
          strm << vec_eigenvalues.at(i).at(j);
        }
      }
      for (std::size_t i = 0u; i < vec_mw_coordinatesMatrix.size(); ++i)
      {
        for (std::size_t j = 0u; j < vec_mw_coordinatesMatrix.at(0u).size(); ++j)
        {
          strm << vec_mw_coordinatesMatrix.at(i).at(j);
        }
      }
      //for (auto const& a : vec_modes)
      //{
      //  for (auto const& b : a)
      //  {
      //    strm << b;
      //  }
      //}
      //for (auto const& a : vec_eigenvectors)
      //{
      //  for (auto const& b : a)
      //  {
      //    strm << b;
      //  }
      //}
      //for (auto const& a : vec_eigenvalues)
      //{
      //  for (auto const& b : a)
      //  {
      //    strm << b;
      //  }
      //}
      //for (auto const& a : vec_mw_coordinatesMatrix)
      //{
      //  for (auto const& b : a)
      //  {
      //    strm << b;
      //  }
      //}
      //strm << modes << eigenvectors << eigenvalues << mw_coordinatesMatrix;

      return strm;
    }

    //**overload for >> operator*/
    template<class Strm>
    friend scon::binary_stream<Strm>& operator>> (scon::binary_stream<Strm>& strm, pca::PrincipalComponentRepresentation& pca)
    {
      std::array<std::size_t, 8u> sizes;
      // sizes
      for (auto& s : sizes)
        strm >> s;

      std::vector<std::vector<float_type>> modes, eigenvec, eigenval, trajectory_mw;
      modes = std::vector<std::vector<float_type>>(sizes[0],std::vector<float_type>(sizes[1], -1.));
      eigenvec = std::vector<std::vector<float_type>>(sizes[2], std::vector<float_type>(sizes[3], -1.));
      eigenval = std::vector<std::vector<float_type>>(sizes[4], std::vector<float_type>(sizes[5], -1.));
      trajectory_mw = std::vector<std::vector<float_type>>(sizes[6], std::vector<float_type>(sizes[7], -1.));
      for (std::size_t i = 0u; i < modes.size(); ++i)
      {
        for (std::size_t j = 0u; j < modes.at(0u).size(); ++j)
        {
          strm >> modes.at(i).at(j);
        }
      }
      for (std::size_t i = 0u; i < eigenvec.size(); ++i)
      {
        for (std::size_t j = 0u; j < eigenvec.at(0u).size(); ++j)
        {
          strm >> eigenvec.at(i).at(j);
        }
      }
      for (std::size_t i = 0u; i < eigenval.size(); ++i)
      {
        for (std::size_t j = 0u; j < eigenval.at(0u).size(); ++j)
        {
          strm >> eigenval.at(i).at(j);
        }
      }
      for (std::size_t i = 0u; i < trajectory_mw.size(); ++i)
      {
        for (std::size_t j = 0u; j < trajectory_mw.at(0u).size(); ++j)
        {
          strm >> trajectory_mw.at(i).at(j);
        }
      }
      //for (auto const& a : modes)
      //{
      //  for (auto const& b : a)
      //  {
      //    strm >> b;
      //  }
      //}
      //for (auto const& a : eigenvec)
      //{
      //  for (auto const& b : a)
      //  {
      //    strm >> b;
      //  }
      //}
      //for (auto const& a : eigenval)
      //{
      //  for (auto const& b : a)
      //  {
      //    strm >> b;
      //  }
      //}
      //for (auto const& a : trajectory_mw)
      //{
      //  for (auto const& b : a)
      //  {
      //    strm >> b;
      //  }
      //}

      //
      //
      pca.modes = Matrix_Class(modes.size(),modes.at(0u).size(),-1.);
      pca.eigenvectors = Matrix_Class(eigenvec.size(), eigenvec.at(0u).size(), -1.);
      pca.eigenvalues = Matrix_Class(eigenval.size(), eigenval.at(0u).size(), -1.);
      pca.mw_coordinatesMatrix = Matrix_Class(trajectory_mw.size(), trajectory_mw.at(0u).size(), -1.);
      for (std::size_t i = 0u; i < pca.modes.rows(); ++i)
      {
        for (std::size_t j = 0u; j < pca.modes.cols(); ++j)
        {
          pca.modes(i,j) = modes.at(i).at(j);
        }
      }
      for (std::size_t i = 0u; i < pca.eigenvectors.rows(); ++i)
      {
        for (std::size_t j = 0u; j < pca.eigenvectors.cols(); ++j)
        {
          pca.eigenvectors(i, j) = eigenvec.at(i).at(j);
        }
      }
      for (std::size_t i = 0u; i < pca.eigenvalues.rows(); ++i)
      {
        for (std::size_t j = 0u; j < pca.eigenvalues.cols(); ++j)
        {
          pca.eigenvalues(i, j) = eigenval.at(i).at(j);
        }
      }
      for (std::size_t i = 0u; i < pca.mw_coordinatesMatrix.rows(); ++i)
      {
        for (std::size_t j = 0u; j < pca.mw_coordinatesMatrix.cols(); ++j)
        {
          pca.mw_coordinatesMatrix(i, j) = trajectory_mw.at(i).at(j);
        }
      }
      // DEBUG
      // std::cout << "Modes:\n";
      // std::cout << pca.getModes() << "\n";
      // std::cout << "Eigenval:\n" << "\n";
      // std::cout << pca.getEigenvalues() << "\n";
      // std::cout << "Eigenvec:\n" << "\n";
      // std::cout << pca.getEigenvectors() << "\n";
      // std::cout << "MWTrRaj:\n" << "\n";
      // std::cout << pca.getMWTrajectoryMatrix() << "\n";
      // std::cout << "#####\n" << "\n";
      return strm;
    }
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


