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

  Matrix_Class unmassweightedStdDevFromMWPCAeigenvalues(Matrix_Class const& massVector, Matrix_Class const& pcaEigenvalues, Matrix_Class const& pcaEigenvectors, std::vector<size_t> const& subDims)
  {
    Matrix_Class assocRedMasses = calculateReducedMassOfPCAModes(massVector, pcaEigenvalues, pcaEigenvectors, subDims);
    std::cout << "SDEBUG assocRedMasses: " << assocRedMasses << std::endl;
    if (assocRedMasses.cols() != 1u || pcaEigenvalues.cols() != 1u || pcaEigenvalues.rows() != assocRedMasses.rows())
    {
      throw std::logic_error("Cannot un-massweight PCA Modes: Dimensionality of matrices is wrong. Aborting!");
      return Matrix_Class();
    }
    Matrix_Class unweightedPcaStdDevs(pcaEigenvalues);
    for (std::size_t i = 0u; i < assocRedMasses.rows(); ++i)
    {
      const double toBeDivided = pcaEigenvalues(i, 0u);
      const double divisor = assocRedMasses(i, 0u);
      unweightedPcaStdDevs(i, 0u) = std::sqrt(toBeDivided) / std::sqrt(divisor);
      std::cout << "SDEBUG: sqrt(" << toBeDivided << ")/sqrt(" << divisor << ")= " << unweightedPcaStdDevs(i, 0u) << std::endl;
    }
    std::cout << "Debug unweighted std-devs:\n" << unweightedPcaStdDevs << "\n";
    return unweightedPcaStdDevs;
  }

  Matrix_Class calculateReducedMassOfPCAModes(Matrix_Class const& massVector, Matrix_Class const& pca_eigenvalues, Matrix_Class const& pca_eigenvectors, std::vector<size_t> const& subDims)
  {
    Matrix_Class assocRedMasses(pca_eigenvalues.rows(), 1u);
    for (std::size_t i = 0; i < pca_eigenvalues.rows(); i++)
    {
      if (subDims == std::vector<size_t>() || std::find(subDims.begin(), subDims.end(), i) != subDims.end())
      {
        if (massVector.cols() == 1u && massVector.rows() == pca_eigenvalues.rows())
        {
          std::cout << "....................\n";
          std::cout << "Debug: Mode " << i << std::endl;
          std::cout << "Debug: eigenvalues " << pca_eigenvalues << std::endl;
          // Assoc red mass of each mode via https://physics.stackexchange.com/questions/401370/normal-modes-how-to-get-reduced-masses-from-displacement-vectors-atomic-masses
          double A___normalizationThisEigenvector = 0.;
          double inv_red_mass = 0.0;
          for (std::size_t j = 0; j < pca_eigenvectors.cols(); j++)
          {
            //Each column is one eigenvector
            const double squaredEigenvecValue = pca_eigenvectors(i, j) * pca_eigenvectors(i, j);
            A___normalizationThisEigenvector += squaredEigenvecValue;
            std::cout << "Debug: A___normalizationThisEigenvector " << A___normalizationThisEigenvector << std::endl;
            const double currentMass = massVector(i, 0u);
            std::cout << "Debug: currentMass " << currentMass << std::endl;
            inv_red_mass += A___normalizationThisEigenvector / currentMass;
            std::cout << "Debug: inv_red_mass currently  " << inv_red_mass << std::endl;
          }
          //
          const double red_mass = 1.0 / inv_red_mass;
          std::cout << "Debug: red_mass " << red_mass << std::endl;
          assocRedMasses(i, 0u) = red_mass;
        }
      }
    }
    return assocRedMasses;
  }


}
