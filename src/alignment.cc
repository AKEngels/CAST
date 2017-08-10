#include "alignment.h"

namespace align
{
  using namespace matop;
  float_type drmsd_calc(coords::Coordinates const& input, coords::Coordinates const& ref)
  {
    if (input.atoms().size() != ref.atoms().size()) throw std::logic_error("Number of atoms of structures passed to drmsd_calc to not match.");

    float_type value = 0;
    for (size_t i = 0; i < input.atoms().size(); i++)
    {
      for (size_t j = 0; j < i; j++)
      {
        float_type holder = sqrt((ref.xyz(i).x() - ref.xyz(j).x()) * (ref.xyz(i).x() - ref.xyz(j).x()) + (ref.xyz(i).y() - ref.xyz(j).y())* (ref.xyz(i).y() - ref.xyz(j).y()) + (ref.xyz(i).z() - ref.xyz(j).z()) * (ref.xyz(i).z() - ref.xyz(j).z()));
        float_type holder2 = sqrt((input.xyz(i).x() - input.xyz(j).x()) * (input.xyz(i).x() - input.xyz(j).x()) + (input.xyz(i).y() - input.xyz(j).y())* (input.xyz(i).y() - input.xyz(j).y()) + (input.xyz(i).z() - input.xyz(j).z()) * (input.xyz(i).z() - input.xyz(j).z()));
        value += (holder2 - holder) * (holder2 - holder);
      }
    }
    return sqrt(value / (double)(input.atoms().size() * (input.atoms().size() + 1u)));
  }

  float_type holmsander_calc(coords::Coordinates const& input, coords::Coordinates const& ref, double holmAndSanderDistance)
  {
    if (input.atoms().size() != ref.atoms().size()) throw std::logic_error("Number of atoms of structures passed to drmsd_calc to not match.");

    float_type value(0);
    for (size_t i = 0; i < input.atoms().size(); i++) {
      for (size_t j = 0; j < i; j++)
      {
        float_type holder = sqrt((ref.xyz(i).x() - ref.xyz(j).x()) * (ref.xyz(i).x() - ref.xyz(j).x()) + (ref.xyz(i).y() - ref.xyz(j).y())* (ref.xyz(i).y() - ref.xyz(j).y()) + (ref.xyz(i).z() - ref.xyz(j).z()) * (ref.xyz(i).z() - ref.xyz(j).z()));
        float_type holder2 = sqrt((input.xyz(i).x() - input.xyz(j).x()) * (input.xyz(i).x() - input.xyz(j).x()) + (input.xyz(i).y() - input.xyz(j).y())* (input.xyz(i).y() - input.xyz(j).y()) + (input.xyz(i).z() - input.xyz(j).z()) * (input.xyz(i).z() - input.xyz(j).z()));
        value += abs(holder2 - holder) * exp(-1 * (holder2 + holder)*(holder2 + holder) / (4 * holmAndSanderDistance * holmAndSanderDistance)) / (holder2 + holder);
      }
    }
    return value;
  }

  coords::Coordinates kabschAligned(coords::Coordinates const& inputCoords, coords::Coordinates const& reference, bool centerOfMassAlign)
  {
    coords::Coordinates output(inputCoords);
    if (centerOfMassAlign) centerOfMassAlignment(output);
    kabschAlignment(output, reference);
    return output;
  }

  void kabschAlignment(coords::Coordinates& inputCoords, coords::Coordinates const& reference, bool centerOfMassAlign)
  {
    if (centerOfMassAlign)
    {
      centerOfMassAlignment(inputCoords);
    }

    Matrix_Class input = transfer_to_matr(inputCoords);
    Matrix_Class ref = transfer_to_matr(reference);

    Matrix_Class c(input * transposed(ref));
    //Creates Covariance Matrix

    Matrix_Class s, V, U;
    c.singular_value_decomposition(U, s, V);

    Matrix_Class unit = Matrix_Class::Identity(c.rows(), c.rows());
    if ((c.det_sign() < 0)) //Making sure that U will do a proper rotation (rows/columns have to be right handed system)
    {
      unit(2, 2) = -1;
    }
    transpose(U);
    unit = unit * U;
    unit = V * unit;
    input = unit * input;

    inputCoords.set_xyz(transfer_to_3DRepressentation(input));
  }

  void centerOfMassAlignment(coords::Coordinates& coords_in)
  {
    coords::Cartesian_Point com_ref = coords_in.center_of_mass();
    coords_in.move_all_by(-com_ref, true);
  }
}


void alignment(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords)
{
  using namespace matop;
  using namespace align;

  coords::Coordinates coordsReferenceStructure(coords), coordsTemporaryStructure(coords);

  // Check if reference structure is in range
  if (Config::get().alignment.reference_frame_num >= ci->size()) throw std::runtime_error("Reference frame number in ALIGN task is bigger than number of frames in the input structure ensemble.");

  auto temporaryPESpoint = ci->PES()[Config::get().alignment.reference_frame_num].structure.cartesian;

  //Alignment to external reference frame (different file)
  if (!Config::get().alignment.align_external_file.empty())
  {
    std::unique_ptr<coords::input::format> externalReferenceStructurePtr(coords::input::new_format());
    coords::Coordinates externalReferenceStructure(externalReferenceStructurePtr->read(Config::get().alignment.align_external_file));
    if (Config::get().alignment.reference_frame_num >= externalReferenceStructurePtr->PES().size())
    {
      throw std::out_of_range("Requested reference frame number not within reference structure ensemble.");
    }
    temporaryPESpoint = externalReferenceStructurePtr->PES()[Config::get().alignment.reference_frame_num].structure.cartesian;
  }
  //Constructs two coordinate objects and sets reference frame according to INPUTFILE
  coordsReferenceStructure.set_xyz(temporaryPESpoint);

  //Construct and Allocate arrays for stringoutput (necessary for OpenMP)
  double mean_value = 0;
  std::string *hold_str, *hold_coords_str;
  hold_str = new std::string[ci->size()];
  hold_coords_str = new std::string[ci->size()];

  //Perform translational alignment for reference frame
  if (Config::get().alignment.traj_align_translational)
  {
    centerOfMassAlignment(coordsReferenceStructure);
  }

  // Output text
  if (Config::get().general.verbosity > 2U) std::cout << "ALIGN preparations done. Starting actual alignment.\n";

#ifdef _OPENMP
  if (Config::get().general.verbosity > 3U) std::cout << "Using openMP for alignment.\n";
  auto const n_omp = static_cast<std::ptrdiff_t>(ci->size());
#pragma omp parallel for firstprivate(coordsReferenceStructure, coordsTemporaryStructure) reduction(+:mean_value) shared(hold_coords_str, hold_str)
  for (std::ptrdiff_t i = 0; i < n_omp; ++i)
#else
  for (std::size_t i = 0; i < ci->size(); ++i)
#endif
  {
    if (i != static_cast<std::ptrdiff_t>(Config::get().alignment.reference_frame_num))
    {
      auto temporaryPESpoint2 = ci->PES()[i].structure.cartesian;
      coordsTemporaryStructure.set_xyz(temporaryPESpoint2);
      //Create temporary objects for current frame

      if (Config::get().alignment.traj_align_translational)
      {
        centerOfMassAlignment(coordsTemporaryStructure);
      }
      if (Config::get().alignment.traj_align_rotational)
      {
        kabschAlignment(coordsTemporaryStructure, coordsReferenceStructure, false);
      }

      if (Config::get().alignment.traj_print_bool)
      {
        if (Config::get().alignment.dist_unit == 0)
          //RMSD
        {
          std::stringstream temporaryStringstream;
          double currentRootMeanSquareDevaition = root_mean_square_deviation(coordsTemporaryStructure.xyz(), coordsReferenceStructure.xyz());
          temporaryStringstream << std::setw(13) << i << " ";
          temporaryStringstream << std::setw(13) << currentRootMeanSquareDevaition << "\n";
          mean_value += currentRootMeanSquareDevaition;
          hold_str[i] = temporaryStringstream.str();
        }
        else if (Config::get().alignment.dist_unit == 1)
          //dRMSD
        {
          std::stringstream temporaryStringstream;
          temporaryStringstream << i << " ";
          double value = (double)drmsd_calc(coordsTemporaryStructure, coordsReferenceStructure);
          temporaryStringstream << std::setw(13) << value << "\n";
          mean_value += value;
          hold_str[i] = temporaryStringstream.str();
        }
        else if (Config::get().alignment.dist_unit == 2)
          //Holm&Sander Distance
        {
          std::stringstream temporaryStringstream;
          double value = (double)holmsander_calc(coordsTemporaryStructure, coordsReferenceStructure, Config::get().alignment.holm_sand_r0);
          temporaryStringstream << std::setw(13) << i << " " << value << "\n";
          mean_value += value;
          hold_str[i] = temporaryStringstream.str();
        }
      }
      //Molecular distance measure calculation

      std::stringstream hold_coords;
      hold_coords << coordsTemporaryStructure;
      hold_coords_str[i] = hold_coords.str();
      //Formatted string-output
    }

    else if (i == static_cast<std::ptrdiff_t>(Config::get().alignment.reference_frame_num))
    {
      std::stringstream hold_coords;
      hold_coords << coordsReferenceStructure;
      hold_coords_str[i] = hold_coords.str();
      //Formatted string-output (first to array because of OpenMP parallelization)
    }
  }

  std::ofstream distance(coords::output::filename("_distances").c_str(), std::ios::app);
  std::ofstream outputstream(coords::output::filename("_aligned").c_str(), std::ios::app);

  if (Config::get().general.verbosity > 2U) std::cout << "Alignment done. Writing structures to file.\n";

  for (size_t i = 0; i < ci->size(); i++)
  {
    if (Config::get().alignment.traj_print_bool)
    {
      distance << hold_str[i];
    }
    outputstream << hold_coords_str[i];
  }
  distance << "\n";
  distance << "Mean value: " << (mean_value / (double)(ci->size() - 1)) << "\n";
  //Formatted string-output

  delete[] hold_str;
  delete[] hold_coords_str;
  //Cleaning Up
}
