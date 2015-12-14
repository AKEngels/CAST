#pragma once 

#if defined _OPENMP
  #include <omp.h>
#endif
#include <memory>
#include <string>
#include <iostream>
#include <sstream>
#include "scon_utility.h"
#include "coords.h"

namespace coords
{
  namespace input
  {
    // format types
    struct types { enum T { TINKER, PDB }; };

    inline types::T get_type (std::string const &str)
    {
      if (str.find("PDB") != str.npos) return types::PDB;
      else return types::TINKER;
    }
    // new format creator
    input::format* new_format (void);

    class format
    {
    protected:
      format (void) { }
      coords::Ensemble_PES input_ensemble;
    public:
      virtual ~format (void) { }
      // Number of Atoms (N) and Structures (M)
      std::size_t atoms (void) const { return input_ensemble.size() > 0 ? input_ensemble.back().size() : 0U; }
      std::size_t size (void) const { return input_ensemble.size(); }
      // Read Structure and return coordinates
      virtual Coordinates read (std::string) = 0;
      // Get structure i
      PES_Point structure (std::size_t const i=0u) const { return input_ensemble[i]; }
      Ensemble_PES const & PES (void) const { return input_ensemble; }
      Ensemble_PES & PES(void) { return input_ensemble; }
      // Iterators
      Ensemble_PES::iterator begin (void) { return input_ensemble.begin(); }
      Ensemble_PES::iterator end (void) { return input_ensemble.end(); }
      Ensemble_PES::const_iterator begin (void) const { return input_ensemble.begin(); }
      Ensemble_PES::const_iterator end (void) const { return input_ensemble.end(); }
    };

    namespace formats
    {
      class tinker
        : public coords::input::format
      {
      public:
        Coordinates read (std::string);
      private:
        struct line
        {
          Cartesian_Point position;
          std::vector<std::size_t> bonds;
          std::size_t index, tinker_type;
          std::string symbol;
          line (void) 
            : bonds(7U,0U), index(0U),
            tinker_type(0U) {}
          friend std::istream & operator>> (std::istream &is, line &l)
          {
            is >> l.index >> l.symbol >> l.position.x() >> l.position.y() >> l.position.z() >> l.tinker_type 
              >> l.bonds[0u] >> l.bonds[1u] >> l.bonds[2u] >> l.bonds[3u] >> l.bonds[4u] >> l.bonds[5u] >> l.bonds[6u];
            return is;
          }
        };
        
      };
    }

  }

  namespace output
  {




    class format
    {
    protected:
      Coordinates const & ref;
      format& operator= (format const &);
      format (Coordinates const &coord_obj) : ref(coord_obj) {}
    public:
      virtual void to_stream (std::ostream&) const = 0;
    };

    inline std::ostream& operator<< (std::ostream &stream, format const & out_format)
    {
      out_format.to_stream(stream);
      return stream;
    }

    namespace formats
    {

      class tinker
        : public output::format
      {
        tinker& operator= (tinker const &);
      public:
        tinker (Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream (std::ostream&) const;
        static std::string filename (std::string postfix)
        {
          return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(".arc")).get_unique_path();
        }
      };

      class xyz
        : public output::format
      {
        xyz& operator= (xyz const &);
      public:
        xyz (Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream (std::ostream&) const;
        static std::string filename (std::string postfix)
        {
          return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(".xyz")).get_unique_path();
        }
      };

      class xyz_mopac
        : public output::format
      {
        xyz_mopac& operator= (xyz_mopac const &);
      public:
        xyz_mopac(Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream(std::ostream&) const;
        static std::string filename(std::string postfix)
        {
          return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(".xyz")).get_unique_path();
        }
      };

	  class xyz_mopac7
		  : public output::format
	  {
		  xyz_mopac7 operator= (xyz_mopac7 const &);
	  public:
		  xyz_mopac7(Coordinates const &coord_obj) : output::format(coord_obj) {}
		  void to_stream(std::ostream&) const;

		  static std::string filename(std::string postfix)
		  {
			  return std::string(Config::get().general.outputFilename).append(postfix).append(".xyz");
		  }
	  };
      class moldenxyz
        : public output::format
      {
        xyz& operator= (xyz const &);
      public:
        moldenxyz(Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream(std::ostream&) const;
        static std::string filename(std::string postfix)
        {
          return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(".xyz")).get_unique_path();
        }
      };

      class zmatrix
        : public output::format
      {
        tinker& operator= (tinker const &);
      public:
        zmatrix (Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream (std::ostream&) const;
        static std::string filename (std::string postfix)
        {
          return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(".zm")).get_unique_path();
        }
      };

    }

    inline std::string filename(std::string postfix = "", std::string extension = "")
    {
      if (extension.length() > 0)
      {
        return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(extension)).get_unique_path();
      }
      if (Config::get().general.output == config::output_types::TINKER) return formats::tinker::filename(postfix);
      else if (Config::get().general.output == config::output_types::XYZ) return formats::xyz::filename(postfix);
      else if (Config::get().general.output == config::output_types::MOLDEN) return formats::moldenxyz::filename(postfix);
      else if (Config::get().general.output == config::output_types::ZMATRIX) return formats::zmatrix::filename(postfix);
      return scon::StringFilePath(Config::get().general.outputFilename).get_unique_path();
    }

  }


}
