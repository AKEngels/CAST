// For more information see matop.h

#include "matop.h"
namespace matop
{
  /////////////////////////////////////
  //                                ///
  //  S P E C I F I C   T A S K S   ///
  //                                ///
  /////////////////////////////////////

  Matrix_Class transfer_to_matr(coords::Coordinates const& in)
  {
    Matrix_Class out_mat(in.size(), 3u);
    for (size_t l = 0; l < in.size(); l++)
    {
      coords::cartesian_type tempcoord2;
      tempcoord2 = in.xyz(l);
      out_mat(l, 0) = tempcoord2.x();
      out_mat(l, 1) = tempcoord2.y();
      out_mat(l, 2) = tempcoord2.z();
    }
    return transposed(out_mat);
  }

  Matrix_Class transfer_to_matr_internal(coords::Coordinates const& in)
  {
    const size_t sizer = in.size();
    Matrix_Class out_mat(sizer, 3u);
    for (size_t l = 0; l < in.size(); l++)
    {
      out_mat(l, 0) = in.intern(l).radius();
      out_mat(l, 1) = in.intern(l).inclination().radians();
      out_mat(l, 2) = in.intern(l).azimuth().radians();
    }
    return transposed(out_mat);
  }

  Matrix_Class transform_coordinates(coords::Coordinates& input)
  {
	  Matrix_Class output(3, input.size());
	  for (size_t l = 0; l < input.size(); l++)
	  {
		  (output)((l), 0) = input.xyz(l).x();
		  (output)((l), 1) = input.xyz(l).y();
		  (output)((l), 2) = input.xyz(l).z();
	  }
    return output;
  }

  coords::Representation_3D transfer_to_3DRepressentation(Matrix_Class const& input)
  {
    coords::Representation_3D tempcoord1;

    for (size_t i = 0; i < input.cols(); i++)
    {
      coords::Cartesian_Point tempcoord2(input(0u, i), input(1u, i), input(2u, i));
      tempcoord1.push_back(tempcoord2);
    }
    return tempcoord1;
  }

  coords::Representation_Internal transfer_to_internalRepressentation(Matrix_Class const& input)
  {
    coords::Representation_Internal tempcoord1;

    for (size_t i = 0; i < input.cols(); i++)
    {
      coords::internal_type tempcoord2(input(0u, i), scon::ang<float_type>(input(1u, i)), scon::ang<float_type>(input(2u, i)));
      tempcoord1.push_back(tempcoord2);
    }
    return tempcoord1;
  }

  void massweight(Matrix_Class& input, coords::Coordinates const& coords, bool to_meter, std::vector<size_t> atomsThatAreUsed)
  //coords are reference coord object to get atomic masses from
  //boolean controls whether coords should also be multiplied with 10e-10 to convert angstrom to meters)
  {
    if (atomsThatAreUsed.empty())
    {
      for (size_t i = 0; i < input.rows(); i = i + 3)
      {
        double temp = sqrt(coords.atoms(i / 3u).mass() * 1.6605402 * 10e-27);
        if (to_meter)
        {
          temp *= 10e-10;
        }
        for (size_t j = 0; j < input.cols(); j++)
        {
          for (size_t k = 0; k < 3; k++)
          {
            (input)(i + k, j) *= temp;
          }
        }
      }
    }
    else
    {
      for (size_t i = 0; i < input.rows(); i = i + 3)
      {
        double temp = sqrt(coords.atoms(atomsThatAreUsed[i / 3u]).mass() * 1.6605402 * 10e-27);
        if (to_meter)
        {
          temp *= 10e-10;
        }
        for (size_t j = 0; j < input.cols(); j++)
        {
          for (size_t k = 0; k < 3; k++)
          {
            (input)(i + k, j) *= temp;
          }
        }
      }
    }
  }

  void undoMassweight(Matrix_Class& input, coords::Coordinates const& coords, bool to_meter, std::vector<size_t> atomsThatAreUsed)
  {
    if (atomsThatAreUsed.empty())
    {
      for (size_t i = 0; i < input.rows(); i = i + 3)
      {
        double temp = sqrt(coords.atoms(i / 3u).mass() * 1.6605402 * 10e-27);
        if (to_meter)
        {
          temp *= 10e-10;
        }
        for (size_t j = 0; j < input.cols(); j++)
        {
          for (size_t k = 0; k < 3; k++)
          {
            (input)(i + k, j) /= temp;
          }
        }
      }
    }
    else
    {
      for (size_t i = 0; i < input.rows(); i = i + 3)
      {
        double temp = sqrt(coords.atoms(atomsThatAreUsed[i / 3u]).mass() * 1.6605402 * 10e-27);
        if (to_meter)
        {
          temp *= 10e-10;
        }
        for (size_t j = 0; j < input.cols(); j++)
        {
          for (size_t k = 0; k < 3; k++)
          {
            (input)(i + k, j) /= temp;
          }
        }
      }
    }
  }

  Matrix_Class transformToOneline(coords::Coordinates const& coords, std::vector<size_t> const& includedAtoms, bool internalCoordinates)
  {
    //First, some range checks
    if (/*if not all atoms*/ includedAtoms.size() != 0 && includedAtoms[includedAtoms.size() - 1u] > coords.atoms().size())
    {
      throw std::runtime_error("Truncation number is greater than the total number of atoms.");
    }
    else if (internalCoordinates && includedAtoms[0] < 3u)
    {
      throw std::runtime_error("Dihedral with index < 3 specified.");
    }

    if (internalCoordinates)
    {
      // Matrix Size
      Matrix_Class transformed_matrix(1u, (includedAtoms.size() * 2));

      size_t j = 0;
      size_t quicksearch_dih = 0;
      for (size_t i = 0; i < coords.atoms().size(); i++)
      {
        size_t keeper = 0;
        bool checker_dih = false;
        for (size_t l = quicksearch_dih; l < Config::get().PCA.pca_internal_dih.size(); l++)
        {
          if (Config::get().PCA.pca_internal_dih.size() != 0)
          {
            if (Config::get().PCA.pca_internal_dih[l] == i)
            {
              checker_dih = true;
              quicksearch_dih++;
              break;
            }
          }
        }
        if (checker_dih)
        {
          transformed_matrix(0, j + keeper) = cos(coords.intern(i).azimuth().radians());
          transformed_matrix(0, j + keeper + 1) = sin(coords.intern(i).azimuth().radians());
          keeper += 2;
        }
        j = j + keeper;
      }
      return transformed_matrix;
    }
    else
    {
      if (includedAtoms.size() != 0u)
      {
        Matrix_Class transformed_matrix;
        transformed_matrix = Matrix_Class(1u, (3 * includedAtoms.size()));

        int j = 0;
        size_t quicksearch = 0;
        for (size_t i = 0; i < coords.atoms().size(); i++) //iterate over atoms
        {
          bool checker = false;
          for (size_t l = quicksearch; l < includedAtoms.size(); l++) //iterate over vector of atoms to account for
          {
            if (includedAtoms[l] == i)
            {
              checker = true;
              quicksearch++;
              break;
            }
          }
          if (checker)
          {
            transformed_matrix(0, j + 0u) = coords.xyz(i).x();
            transformed_matrix(0, j + 1u) = coords.xyz(i).y();
            transformed_matrix(0, j + 2u) = coords.xyz(i).z();
            j = j + 3u;
          }
        }
        return transformed_matrix;
      }
      else
      {
          Matrix_Class transformed_matrix(1, (coords.atoms().size() * 3u));
          int j = 0;
          for (size_t i = 0; i < coords.atoms().size(); i++)
          {
            for (size_t k = 0; k < 3; k++)
            {
              transformed_matrix(0, j + 0u) = coords.xyz(i).x();
              transformed_matrix(0, j + 1u) = coords.xyz(i).y();
              transformed_matrix(0, j + 2u) = coords.xyz(i).z();
            }
            j = j + 3u;
          }
          return transformed_matrix;
      }
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                    W R A P P E R F U N C T I O N S                         //
////////////////////////////////////////////////////////////////////////////////



