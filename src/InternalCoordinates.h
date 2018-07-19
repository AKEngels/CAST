#ifndef H_INTERNAL_COORDINATES
#define H_INTERNAL_COORDINATES

#include<array>

#include "coords.h"
#include "scon_mathmatrix.h"

namespace InternalCoordinates {

  class AbstractGeometryObserver {
  public:
    virtual void update() = 0;
  protected:
    virtual void notify() = 0;
  };

  class AbstractRotatorListener {
  public:
    virtual void setAllFlag() = 0;
  };

  class RotatorObserver : public AbstractGeometryObserver {
  public:
    void setNewRotator(std::shared_ptr<AbstractRotatorListener> const rotator);
    void update() override;
  private:
    void notify() override;
    std::shared_ptr<AbstractRotatorListener> rotator;
  };

  class CartesiansForInternalCoordinates {
  public:
    template<typename T>
    CartesiansForInternalCoordinates(T&& cartesianCoordinates);

    void registerObserver(std::shared_ptr<RotatorObserver> const observer);

    template<typename T>
    void setCartesianCoordnates(T&& newCartesianCoordinates);
    coords::Representation_3D & getCartesianCoordnates() { return cartesianCoordinates; };
    coords::Representation_3D const& getCartesianCoordnates() const { return cartesianCoordinates; };

    coords::r3 & at(std::size_t const index) { return cartesianCoordinates.at(index); }
    coords::r3 const& at(std::size_t const index) const { return cartesianCoordinates.at(index); }

    std::size_t size() const { return cartesianCoordinates.size(); }

  private:
    coords::Representation_3D cartesianCoordinates;
    std::vector<std::shared_ptr<AbstractGeometryObserver>> observerList;
    void notify();
  };

  inline coords::Representation_3D sliceCartesianCoordinates(CartesiansForInternalCoordinates const& cartesians, std::vector<std::size_t> const& indexVector) {
    coords::Representation_3D slicedCoordinates;
    for (auto const& index : indexVector) {
      slicedCoordinates.emplace_back(cartesians.at(index-1u));
    }
    return slicedCoordinates;
  }

  struct InternalCoordinate {
    virtual coords::float_type val(CartesiansForInternalCoordinates const& cartesians) const = 0;
    virtual std::vector<coords::float_type> der_vec(CartesiansForInternalCoordinates const& cartesians) const = 0;
    virtual coords::float_type hessian_guess(CartesiansForInternalCoordinates const& cartesians) const = 0;
    virtual std::string info(CartesiansForInternalCoordinates const & cartesians) const = 0;
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

    coords::float_type val(CartesiansForInternalCoordinates const& cartesians) const override;
    std::pair<coords::r3, coords::r3> der(CartesiansForInternalCoordinates const& cartesians) const;
    std::vector<coords::float_type> der_vec(CartesiansForInternalCoordinates const& cartesians) const override;
    coords::float_type hessian_guess(CartesiansForInternalCoordinates const& cartesians) const override;
    std::string info(CartesiansForInternalCoordinates const& cartesians) const override;
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

    coords::float_type val(CartesiansForInternalCoordinates const& cartesians) const override;
    std::tuple<coords::r3, coords::r3, coords::r3> der(CartesiansForInternalCoordinates const& cartesians) const;
    std::vector<coords::float_type> der_vec(CartesiansForInternalCoordinates const& cartesians) const override;
    coords::float_type hessian_guess(CartesiansForInternalCoordinates const& cartesians) const override;
    std::string info(CartesiansForInternalCoordinates const& cartesians) const override;
  };

  struct DihedralAngle : public InternalCoordinates::InternalCoordinate {
    DihedralAngle(const unsigned int& index_a, const unsigned int& index_b,
      const unsigned int& index_c, const unsigned int& index_d)
      : index_a_{ index_a },
      index_b_{ index_b }, index_c_{ index_c }, index_d_{ index_d } {}

    coords::float_type val(CartesiansForInternalCoordinates const& cartesians) const override;
    std::tuple<coords::r3, coords::r3, coords::r3, coords::r3>
      der(CartesiansForInternalCoordinates const& cartesians) const;
    std::vector<coords::float_type> der_vec(CartesiansForInternalCoordinates const& cartesians) const override;
    coords::float_type hessian_guess(CartesiansForInternalCoordinates const& cartesians) const override;
    std::string info(CartesiansForInternalCoordinates const& cartesians) const override;

    std::size_t index_a_;
    std::size_t index_b_;
    std::size_t index_c_;
    std::size_t index_d_;
  };

  struct OutOfPlane : public DihedralAngle {
    using DihedralAngle::DihedralAngle;

    coords::float_type hessian_guess(CartesiansForInternalCoordinates const& cartesians) const override;
    std::string info(CartesiansForInternalCoordinates const& cartesians) const override;
  };

  struct Translations : public InternalCoordinates::InternalCoordinate {
    Translations(std::vector<std::size_t> const& index_vec)
      : indices_(index_vec) {}

    virtual coords::float_type val(CartesiansForInternalCoordinates const& cartesians) const = 0;
    std::vector<std::size_t> indices_;

    template<typename Func> //TODO remove sys_size!!!!
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

    coords::float_type hessian_guess(CartesiansForInternalCoordinates const& /*cartesians*/) const override {
      return 0.05;
    }

  };

  struct TranslationX : Translations {
    TranslationX(const std::vector<std::size_t>& index_vec)
      : Translations(index_vec) {}

    coords::float_type val(CartesiansForInternalCoordinates const& cartesians) const override {
      auto coord_sum{ 0.0 };
      for (auto& i : indices_) {
        coord_sum += cartesians.at(i - 1u).x();
      }
      return coord_sum / indices_.size();
    }

    std::vector<coords::float_type> der_vec(CartesiansForInternalCoordinates const& rep)const override;
    std::string info(CartesiansForInternalCoordinates const& cartesians) const override;
  };

  struct TranslationY : Translations {
    TranslationY(const std::vector<std::size_t>& index_vec)
      : Translations(index_vec) {}

    coords::float_type val(CartesiansForInternalCoordinates const& cartesians) const override {
      auto coord_sum{ 0.0 };
      for (auto& i : indices_) {
        coord_sum += cartesians.at(i - 1u).y();
      }
      return coord_sum / indices_.size();
    }

    std::vector<coords::float_type> der_vec(CartesiansForInternalCoordinates const& rep)const override;
    std::string info(CartesiansForInternalCoordinates const& cartesians) const override;
  };

  struct TranslationZ : Translations {
    TranslationZ(const std::vector<std::size_t>& index_vec)
      : Translations(index_vec) {}

    coords::float_type val(CartesiansForInternalCoordinates const& cartesians) const override {
      auto coord_sum{ 0.0 };
      for (auto& i : indices_) {
        coord_sum += cartesians.at(i - 1u).z();
      }
      return coord_sum / indices_.size();
    }

    std::vector<coords::float_type> der_vec(CartesiansForInternalCoordinates const& rep)const override;
    std::string info(CartesiansForInternalCoordinates const& cartesians) const override;
  };

  template<typename T>
  inline CartesiansForInternalCoordinates::CartesiansForInternalCoordinates(T && cartesians) : cartesianCoordinates(std::forward<T>(cartesians)) {}

  template<typename T>
  inline void CartesiansForInternalCoordinates::setCartesianCoordnates(T && newCartesianCoordinates) {
    cartesianCoordinates = std::forward<T>(newCartesianCoordinates);
    notify();
  }

  struct Rotator;
  struct RotationA;
  struct RotationB;
  struct RotationC;

  struct Rotations {
    std::shared_ptr<Rotator> rotator;
    std::unique_ptr<InternalCoordinate> rotationA;
  };


  class Rotator : public AbstractRotatorListener, public std::enable_shared_from_this<Rotator> {
  public:

    static std::shared_ptr<Rotator> buildRotator(InternalCoordinates::CartesiansForInternalCoordinates & cartesians, std::vector<std::size_t> const& indexVector){
      auto newInstance = std::make_shared<Rotator>(sliceCartesianCoordinates(cartesians, indexVector), indexVector);
      newInstance->registerCartesians(cartesians);
      return newInstance;
    }

    void setAllFlag()override { updateStoredValues = updateStoredDerivatives = true; }
    bool areValuesUpToDate() { return !updateStoredValues; }
    bool areDerivativesUpToDate() { return !updateStoredDerivatives; }

    std::array<coords::float_type, 3u> const& valueOfInternalCoordinate(const coords::Representation_3D&);
    scon::mathmatrix<coords::float_type> const& rot_der_mat(const coords::Representation_3D&);

    Rotations makeRotations();

    //TODO make it private again
    Rotator(coords::Representation_3D const& reference, std::vector<std::size_t> const& index_vec) :
      reference_{ reference }, indices_{ index_vec },
      rad_gyr_{ radiusOfGyration(reference_) }, updateStoredValues{ true }, updateStoredDerivatives{ true } {}
  private:

    std::vector<scon::mathmatrix<coords::float_type>> rot_der(const coords::Representation_3D&) const;
    coords::float_type radiusOfGyration(const coords::Representation_3D&);


    std::array<coords::float_type, 3u> storedValuesForRotations;
    scon::mathmatrix<coords::float_type> storedDerivativesForRotations;
    bool updateStoredValues;
    bool updateStoredDerivatives;

    void registerCartesians(InternalCoordinates::CartesiansForInternalCoordinates & cartesianCoordinates);

    coords::Representation_3D const reference_;
    std::vector<std::size_t> indices_;
    coords::float_type rad_gyr_;
  };


  struct RotationA : public InternalCoordinate {
    RotationA(std::shared_ptr<Rotator> const rotator) : rotator{ rotator } {}
    virtual coords::float_type val(CartesiansForInternalCoordinates const& cartesians) const override {
      auto const& returnValues = rotator->valueOfInternalCoordinate(cartesians.getCartesianCoordnates());
      return returnValues.at(0);
    }
    virtual std::vector<coords::float_type> der_vec(CartesiansForInternalCoordinates const& cartesians) const override {
      auto const& derivativeMatrix = rotator->rot_der_mat(cartesians.getCartesianCoordnates());
      return derivativeMatrix.col_to_std_vector(0);
    }
    virtual coords::float_type hessian_guess(CartesiansForInternalCoordinates const& cartesians) const override {
      return 0.0;
    }
    virtual std::string info(CartesiansForInternalCoordinates const & cartesians) const override {
      return "";
    }
    std::shared_ptr<Rotator> rotator;
  };


  

  /*class Rotator : public InternalCoordinates::AbstractRotatorListener, public std::enable_shared_from_this<Rotator> {
  public:
    static std::shared_ptr<Rotator> buildInterestedRotator(InternalCoordinates::CartesiansForInternalCoordinates & cartesians);
    void setUpdateFlag()override { updateFlag = true; }
    bool isFlagSet() { return updateFlag; }
  private:
    Rotator() : updateFlag{ false } {}
    void registerCartesians(InternalCoordinates::CartesiansForInternalCoordinates & cartesianCoordinates);
    bool updateFlag;
  };*/
}

#endif