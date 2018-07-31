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

  //Might rather inherit from Representation_3D
  class CartesiansForInternalCoordinates : public coords::Representation_3D {
  public:
    template<typename T>
    CartesiansForInternalCoordinates(T&& cartesianCoordinates); 
    CartesiansForInternalCoordinates() = default;

    void registerObserver(std::shared_ptr<RotatorObserver> const observer);

    template<typename T>
    CartesiansForInternalCoordinates setCartesianCoordnates(T&& newCartesianCoordinates);

  private:
    std::vector<std::shared_ptr<AbstractGeometryObserver>> observerList;
    void notify();
  };

  inline coords::Representation_3D sliceCartesianCoordinates(CartesiansForInternalCoordinates const& cartesians, std::vector<std::size_t> const& indexVector) {
    coords::Representation_3D slicedCoordinates;
    for (auto const& index : indexVector) {
      slicedCoordinates.emplace_back(cartesians.at(index));
    }
    return slicedCoordinates;
  }

  struct InternalCoordinate {
    virtual coords::float_type val(coords::Representation_3D const& cartesians) const = 0;
    virtual std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const = 0;
    virtual coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const = 0;
    virtual std::string info(coords::Representation_3D const & cartesians) const = 0;
  };

  struct BondDistance : public InternalCoordinate {
    template<typename Atom>
    BondDistance(Atom const& atomOne, Atom const& atomTwo)
      : index_a_{ atomOne.atom_serial }, index_b_{ atomTwo.atom_serial },
      elem_a_{ atomOne.element }, elem_b_{ atomTwo.element } {}

    std::size_t index_a_;
    std::size_t index_b_;
    std::string elem_a_;
    std::string elem_b_;

    coords::float_type val(coords::Representation_3D const& cartesians) const override;
    std::pair<coords::r3, coords::r3> der(coords::Representation_3D const& cartesians) const;
    std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override;
    coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
    std::string info(coords::Representation_3D const& cartesians) const override;

    bool operator==(BondDistance const&) const;
  };

  

  struct BondAngle : InternalCoordinate {
    template<typename Atom>
    BondAngle(Atom const& leftAtom, Atom const& middleAtom, Atom const& rightAtom)
      : index_a_{ leftAtom.atom_serial }, index_b_{ middleAtom.atom_serial },
      index_c_{ rightAtom.atom_serial }, elem_a_{ leftAtom.element }, elem_b_{ middleAtom.element }, elem_c_{
      rightAtom.element
    } {}

    std::size_t index_a_;
    std::size_t index_b_;
    std::size_t index_c_;
    std::string elem_a_;
    std::string elem_b_;
    std::string elem_c_;

    coords::float_type val(coords::Representation_3D const& cartesians) const override;
    std::tuple<coords::r3, coords::r3, coords::r3> der(coords::Representation_3D const& cartesians) const;
    std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override;
    coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
    std::string info(coords::Representation_3D const& cartesians) const override;

    bool operator==(BondAngle const&) const;
  };

  struct DihedralAngle : public InternalCoordinates::InternalCoordinate {
    DihedralAngle(const unsigned int& index_a, const unsigned int& index_b,
      const unsigned int& index_c, const unsigned int& index_d)
      : index_a_{ index_a },
      index_b_{ index_b }, index_c_{ index_c }, index_d_{ index_d } {}

    coords::float_type val(coords::Representation_3D const& cartesians) const override;
    std::tuple<coords::r3, coords::r3, coords::r3, coords::r3>
      der(coords::Representation_3D const& cartesians) const;
    std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override;
    coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
    std::string info(coords::Representation_3D const& cartesians) const override;

    std::size_t index_a_;
    std::size_t index_b_;
    std::size_t index_c_;
    std::size_t index_d_;

    bool operator==(DihedralAngle const&) const;
  };

  struct OutOfPlane : public DihedralAngle {
    using DihedralAngle::DihedralAngle;

    coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
    std::string info(coords::Representation_3D const& cartesians) const override;
  };

  struct Translations : public InternalCoordinates::InternalCoordinate {
    Translations(std::vector<std::size_t> const& index_vec)
      : indices_(index_vec) {}

    //virtual coords::float_type val(coords::Representation_3D const& cartesians) const = 0;
    std::vector<std::size_t> indices_;

    template<typename Func>
    coords::Representation_3D der(std::size_t const & sys_size, Func const & coord) const {
      using rep3D = coords::Representation_3D;
      using cp = coords::Cartesian_Point;

      rep3D result(sys_size, cp(0., 0., 0.));

      auto const & s{ indices_.size() };
      for (auto const & i : indices_) {
        result.at(i) = coord(s);
      }
      return result;
    }

    coords::float_type hessian_guess(coords::Representation_3D const& /*cartesians*/) const override {
      return 0.05;
    }

    bool operator==(Translations const&) const;

  };

  struct TranslationX : Translations {
    TranslationX(const std::vector<std::size_t>& index_vec)
      : Translations(index_vec) {}

    coords::float_type val(coords::Representation_3D const& cartesians) const override {
      auto coord_sum{ 0.0 };
      for (auto& i : indices_) {
        coord_sum += cartesians.at(i).x();
      }
      return coord_sum / indices_.size();
    }

    std::vector<coords::float_type> der_vec(coords::Representation_3D const& rep)const override;
    std::string info(coords::Representation_3D const& cartesians) const override;
  };

  struct TranslationY : Translations {
    TranslationY(const std::vector<std::size_t>& index_vec)
      : Translations(index_vec) {}

    coords::float_type val(coords::Representation_3D const& cartesians) const override {
      auto coord_sum{ 0.0 };
      for (auto& i : indices_) {
        coord_sum += cartesians.at(i).y();
      }
      return coord_sum / indices_.size();
    }

    std::vector<coords::float_type> der_vec(coords::Representation_3D const& rep)const override;
    std::string info(coords::Representation_3D const& cartesians) const override;
  };

  struct TranslationZ : Translations {
    TranslationZ(const std::vector<std::size_t>& index_vec)
      : Translations(index_vec) {}

    coords::float_type val(coords::Representation_3D const& cartesians) const override {
      auto coord_sum{ 0.0 };
      for (auto& i : indices_) {
        coord_sum += cartesians.at(i).z();
      }
      return coord_sum / indices_.size();
    }

    std::vector<coords::float_type> der_vec(coords::Representation_3D const& rep)const override;
    std::string info(coords::Representation_3D const& cartesians) const override;
  };

  template<typename T>
  inline CartesiansForInternalCoordinates::CartesiansForInternalCoordinates(T && cartesians) : coords::Representation_3D(std::forward<T>(cartesians)) {}

  template<typename T>
  inline CartesiansForInternalCoordinates CartesiansForInternalCoordinates::setCartesianCoordnates(T && newCartesianCoordinates) {
    notify();
    return *this = std::forward<T>(newCartesianCoordinates);
  }

  struct Rotator;
  struct RotationA;
  struct RotationB;
  struct RotationC;

  struct Rotations {
    std::shared_ptr<Rotator> rotator;
    std::unique_ptr<InternalCoordinate> rotationA;
    std::unique_ptr<InternalCoordinate> rotationB;
    std::unique_ptr<InternalCoordinate> rotationC;
  };


  class Rotator : public AbstractRotatorListener, public std::enable_shared_from_this<Rotator> {
  public:

    static std::shared_ptr<Rotator> buildRotator(InternalCoordinates::CartesiansForInternalCoordinates & cartesians, std::vector<std::size_t> const& indexVector){
      //The ctor should keep being private. Thus make_shared connot be used. If someone has a better idea, go ahead and refactor :)
      auto newInstance = std::shared_ptr<Rotator>(new Rotator(sliceCartesianCoordinates(cartesians, indexVector), indexVector));
      newInstance->registerCartesians(cartesians);
      return newInstance;
    }

    void setAllFlag()override { updateStoredValues = updateStoredDerivatives = true; }
    void requestNewValueEvaluation() { updateStoredValues = true; }
    bool areValuesUpToDate() { return !updateStoredValues; }
    bool areDerivativesUpToDate() { return !updateStoredDerivatives; }

    std::array<coords::float_type, 3u> const& valueOfInternalCoordinate(const coords::Representation_3D&);
    scon::mathmatrix<coords::float_type> const& rot_der_mat(coords::Representation_3D const&);

    Rotations makeRotations();

    bool operator==(Rotator const& other) const;
    
  private:
    Rotator(coords::Representation_3D const& reference, std::vector<std::size_t> const& index_vec) :
      reference_{ reference }, indices_{ index_vec },
      rad_gyr_{ radiusOfGyration(reference_) }, updateStoredValues{ true }, updateStoredDerivatives{ true } {}

    std::vector<scon::mathmatrix<coords::float_type>> rot_der(coords::Representation_3D const&) const;
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
    virtual coords::float_type val(coords::Representation_3D const& cartesians) const override {
      auto const& returnValues = rotator->valueOfInternalCoordinate(cartesians);
      return returnValues.at(0);
    }
    virtual std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override {
      auto const& derivativeMatrix = rotator->rot_der_mat(cartesians);
      return derivativeMatrix.col_to_std_vector(0);
    }
    virtual coords::float_type hessian_guess(coords::Representation_3D const& /*cartesians*/) const override {
      return 0.05;
    }
    virtual std::string info(coords::Representation_3D const & cartesians) const override {
      std::ostringstream oss;
      oss << "Rotation A: " << val(cartesians);
      return oss.str();
    }
    std::shared_ptr<Rotator> rotator;
    bool operator==(RotationA const& other) const {
      return *rotator.get() == *other.rotator.get();
    }
  };


  struct RotationB : public InternalCoordinate {
    RotationB(std::shared_ptr<Rotator> const rotator) : rotator{ rotator } {}
    virtual coords::float_type val(coords::Representation_3D const& cartesians) const override {
      auto const& returnValues = rotator->valueOfInternalCoordinate(cartesians);
      return returnValues.at(1);
    }
    virtual std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override {
      auto const& derivativeMatrix = rotator->rot_der_mat(cartesians);
      return derivativeMatrix.col_to_std_vector(1);
    }
    virtual coords::float_type hessian_guess(coords::Representation_3D const& /*cartesians*/) const override {
      return 0.05;
    }
    virtual std::string info(coords::Representation_3D const & cartesians) const override {
        std::ostringstream oss;
        oss << "Rotation B: " << val(cartesians);
        return oss.str();
    }
    std::shared_ptr<Rotator> rotator;
    bool operator==(RotationB const& other) const {
      return *rotator.get() == *other.rotator.get();
    }
  };

  struct RotationC : public InternalCoordinate {
    RotationC(std::shared_ptr<Rotator> const rotator) : rotator{ rotator } {}
    virtual coords::float_type val(coords::Representation_3D const& cartesians) const override {
      auto const& returnValues = rotator->valueOfInternalCoordinate(cartesians);
      return returnValues.at(2);
    }
    virtual std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override {
      auto const& derivativeMatrix = rotator->rot_der_mat(cartesians);
      return derivativeMatrix.col_to_std_vector(2);
    }
    virtual coords::float_type hessian_guess(coords::Representation_3D const& /*cartesians*/) const override {
      return 0.05;
    }
    virtual std::string info(coords::Representation_3D const & cartesians) const override {
      std::ostringstream oss;
      oss << "Rotation C: " << val(cartesians);
      return oss.str();
    }
    std::shared_ptr<Rotator> rotator;
    bool operator==(RotationC const& other) const {
      return *rotator.get() == *other.rotator.get();
    }
  };
  
}

#endif