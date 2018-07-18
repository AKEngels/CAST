#ifndef H_INTERNAL_COORDINATES
#define H_INTERNAL_COORDINATES

#include "coords.h"

namespace InternalCoordinates {

  struct InternalCoordinate {
    virtual coords::float_type val(coords::Representation_3D const& xyz) const = 0;
    virtual std::vector<coords::float_type> der_vec(coords::Representation_3D const& xyz) const = 0;
    virtual coords::float_type hessian_guess(coords::Representation_3D const& xyz) const = 0;
    virtual std::string info(coords::Representation_3D const & xyz) const = 0;
  };

  struct BondDistance : public InternalCoordinate {
    BondDistance(const unsigned int& index_a, const unsigned int& index_b,
      const std::string& elem_a, const std::string& elem_b)
      : index_a_{ index_a }, index_b_{ index_b },
      elem_a_{ elem_a }, elem_b_{ elem_b } {}

    std::size_t index_a_;
    std::size_t index_b_;
    std::string elem_a_;
    std::string elem_b_;

    coords::float_type val(coords::Representation_3D const& xyz) const override;
    std::pair<coords::r3, coords::r3> der(coords::Representation_3D const& xyz) const;
    std::vector<coords::float_type> der_vec(coords::Representation_3D const& xyz) const override;
    coords::float_type hessian_guess(coords::Representation_3D const& xyz) const override;
    std::string info(coords::Representation_3D const& xyz) const override;
  };


  struct BondAngle : InternalCoordinate {
    BondAngle(const unsigned int& index_a,
      const unsigned int& index_b, const unsigned int& index_c,
      const std::string& elem_a, const std::string& elem_b,
      const std::string& elem_c)
      : index_a_{ index_a }, index_b_{ index_b },
      index_c_{ index_c }, elem_a_{ elem_a }, elem_b_{ elem_b }, elem_c_{
      elem_c
    } {}

    std::size_t index_a_;
    std::size_t index_b_;
    std::size_t index_c_;
    std::string elem_a_;
    std::string elem_b_;
    std::string elem_c_;

    coords::float_type val(coords::Representation_3D const& xyz) const override;
    std::tuple<coords::r3, coords::r3, coords::r3> der(coords::Representation_3D const& xyz) const;
    std::vector<coords::float_type> der_vec(coords::Representation_3D const& xyz) const override;
    coords::float_type hessian_guess(coords::Representation_3D const& xyz) const override;
    std::string info(coords::Representation_3D const& xyz) const override;
  };

  struct DihedralAngle : public InternalCoordinates::InternalCoordinate {
    DihedralAngle(const unsigned int& index_a, const unsigned int& index_b,
      const unsigned int& index_c, const unsigned int& index_d)
      : index_a_{ index_a },
      index_b_{ index_b }, index_c_{ index_c }, index_d_{ index_d } {}

    coords::float_type val(coords::Representation_3D const& xyz) const override;
    std::tuple<coords::r3, coords::r3, coords::r3, coords::r3>
      der(coords::Representation_3D const& xyz) const;
    std::vector<coords::float_type> der_vec(coords::Representation_3D const& xyz) const override;
    coords::float_type hessian_guess(coords::Representation_3D const& xyz) const override;
    std::string info(coords::Representation_3D const& xyz) const override;

    std::size_t index_a_;
    std::size_t index_b_;
    std::size_t index_c_;
    std::size_t index_d_;
  };

  struct OutOfPlane : public DihedralAngle {
    using DihedralAngle::DihedralAngle;

    coords::float_type hessian_guess(coords::Representation_3D const& xyz) const override;
    std::string info(coords::Representation_3D const& xyz) const override;
  };

  struct Translations : public InternalCoordinates::InternalCoordinate {
    Translations(std::vector<std::size_t> const& index_vec)
      : indices_(index_vec) {}

    virtual coords::float_type val(coords::Representation_3D const&) const = 0;
    std::vector<std::size_t> indices_;

    template<typename Func>
    coords::Representation_3D der(std::size_t const & sys_size, Func const & coord) const {
      using rep3D = coords::Representation_3D;
      using cp = coords::Cartesian_Point;

      rep3D result(sys_size, cp(0., 0., 0.));

      auto const & s{ indices_.size() };
      for (auto const & i : indices_) {
        result.at(i - 1) = coord(s);
      }
      return result;
    }

    coords::float_type hessian_guess(coords::Representation_3D const& /*xyz*/) const override {
      return 0.05;
    }

  };

  struct TranslationX : Translations {
    TranslationX(const std::vector<std::size_t>& index_vec)
      : Translations(index_vec) {}

    coords::float_type val(const coords::Representation_3D& xyz) const override {
      auto coord_sum{ 0.0 };
      for (auto& i : indices_) {
        coord_sum += xyz.at(i - 1u).x();
      }
      return coord_sum / indices_.size();
    }

    std::vector<coords::float_type> der_vec(coords::Representation_3D const& rep)const override;
    std::string info(coords::Representation_3D const& xyz) const override;
  };

  struct TranslationY : Translations {
    TranslationY(const std::vector<std::size_t>& index_vec)
      : Translations(index_vec) {}

    coords::float_type val(const coords::Representation_3D& xyz) const override {
      auto coord_sum{ 0.0 };
      for (auto& i : indices_) {
        coord_sum += xyz.at(i - 1u).y();
      }
      return coord_sum / indices_.size();
    }

    std::vector<coords::float_type> der_vec(coords::Representation_3D const& rep)const override;
    std::string info(coords::Representation_3D const& xyz) const override;
  };

  struct TranslationZ : Translations {
    TranslationZ(const std::vector<std::size_t>& index_vec)
      : Translations(index_vec) {}

    coords::float_type val(const coords::Representation_3D& xyz) const override {
      auto coord_sum{ 0.0 };
      for (auto& i : indices_) {
        coord_sum += xyz.at(i - 1u).z();
      }
      return coord_sum / indices_.size();
    }

    std::vector<coords::float_type> der_vec(coords::Representation_3D const& rep)const override;
    std::string info(coords::Representation_3D const& xyz) const override;
  };

}

#endif