#ifndef H_INTERNAL_COORDINATES
#define H_INTERNAL_COORDINATES

#include<array>

#include "coords.h"
#include"ic_atom.h"

namespace scon {
	template <typename T> class mathmatrix;
}

namespace InternalCoordinates {

  class AbstractGeometryObserver {
  public:
    virtual ~AbstractGeometryObserver() = default;
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
    virtual ~RotatorObserver() = default;
    void setNewRotator(std::shared_ptr<AbstractRotatorListener> const rotator);
    void update() override;
  private:
    void notify() override;
    std::shared_ptr<AbstractRotatorListener> rotator;
  };
  
  class temporaryCartesian;

  template<typename CartesianType>
  class CartesiansForInternalCoordinatesImpl : public coords::Container<CartesianType> {
  public:
    template<typename T>
    CartesiansForInternalCoordinatesImpl(T&& cartesianCoordinates); 
    CartesiansForInternalCoordinatesImpl() = default;

    void registerObserver(std::shared_ptr<RotatorObserver> const observer);

    template<typename T>
    CartesiansForInternalCoordinatesImpl& setCartesianCoordnates(T&& newCartesianCoordinates);

    void reset() { notify(); };

  private:
    friend class temporaryCartesian;
    std::vector<std::shared_ptr<AbstractGeometryObserver>> observerList;
    void notify();
  };

  using CartesiansForInternalCoordinates = CartesiansForInternalCoordinatesImpl<coords::Cartesian_Point>;

  class temporaryCartesian{
  public:
    temporaryCartesian(CartesiansForInternalCoordinates & cartesians) : coordinates{ cartesians },
      stolenNotify{ [&cartesians]() { cartesians.notify(); } } {}
    coords::Representation_3D coordinates;
    std::function<void()> stolenNotify;
  };

  inline coords::Representation_3D sliceCartesianCoordinates(CartesiansForInternalCoordinates const& cartesians, std::vector<std::size_t> const& indexVector) {
    coords::Representation_3D slicedCoordinates;
    for (auto const& index : indexVector) {
      slicedCoordinates.emplace_back(cartesians.at(index - 1));
    }
    return slicedCoordinates;
  }

  struct InternalCoordinate {
    virtual coords::float_type val(coords::Representation_3D const& cartesians) const = 0;
    virtual std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const = 0;
    virtual coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const = 0;
    virtual std::string info(coords::Representation_3D const & cartesians) const = 0;
	virtual void makeConstrained() = 0;
	virtual void releaseConstraint() = 0;
    virtual bool is_constrained() const = 0;
    virtual ~InternalCoordinate() = default;
  };

  struct BondDistance : public InternalCoordinate {
    template<typename Atom>
    BondDistance(Atom const& atomOne, Atom const& atomTwo)
      : index_a_{ atomOne.atom_serial - 1u }, index_b_{ atomTwo.atom_serial - 1u },
      elem_a_{ atomOne.element }, elem_b_{ atomTwo.element },
      constrained_{ false /*constraint? *constraint : Config::get().constrained_internals.constrain_bond_lengths*/ }
    {}

    std::size_t index_a_;
    std::size_t index_b_;
    std::string elem_a_;
    std::string elem_b_;

    coords::float_type val(coords::Representation_3D const& cartesians) const override;
    std::pair<coords::r3, coords::r3> der(coords::Representation_3D const& cartesians) const;
    std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override;
    coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
    std::string info(coords::Representation_3D const& cartesians) const override;

	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

    bool operator==(BondDistance const&) const;
    
    bool constrained_;
    virtual bool is_constrained() const override {return constrained_;}

  private:
    bool bothElementsInPeriodOne(ic_atom::period const atomA, ic_atom::period const atomB)const;
    bool oneElementInPeriodOneTheOtherInPeriodTwo(ic_atom::period const atomA, ic_atom::period const atomB)const;
    bool oneElementInPeriodOneTheOtherInPeriodThree(ic_atom::period const atomA, ic_atom::period const atomB)const;
    bool bothElementsInPeriodTwo(ic_atom::period const atomA, ic_atom::period const atomB)const;
    bool oneElementInPeriodTwoTheOtherInPeriodThree(ic_atom::period const atomA, ic_atom::period const atomB)const;
  };

  struct BondAngle : InternalCoordinate {
    template<typename Atom>
    BondAngle(Atom const& leftAtom, Atom const& middleAtom, Atom const& rightAtom)
      : index_a_{ leftAtom.atom_serial - 1u }, index_b_{ middleAtom.atom_serial - 1u },
      index_c_{ rightAtom.atom_serial - 1u }, elem_a_{ leftAtom.element }, elem_b_{ middleAtom.element },
      elem_c_{ rightAtom.element },
      constrained_{ false /*constraint? *constraint : Config::get().constrained_internals.constrain_bond_angles*/ }
    {}

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

	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

    bool operator==(BondAngle const&) const;
    
    bool constrained_;
    virtual bool is_constrained() const override {return constrained_;}
  };

  struct DihedralAngle : public InternalCoordinates::InternalCoordinate {

    template<typename Atom>
    DihedralAngle(Atom const& outerLeftAtom, Atom const& leftAtom,
      Atom const& rightAtom, Atom const& outerRightAtom)
      : index_a_{ outerLeftAtom.atom_serial - 1u },
      index_b_{ leftAtom.atom_serial - 1u }, index_c_{ rightAtom.atom_serial - 1u }, index_d_{ outerRightAtom.atom_serial - 1u },
      constrained_{ false /*constraint? *constraint : Config::get().constrained_internals.constrain_dihedrals*/ }{}
    virtual ~DihedralAngle() = default;

    coords::float_type val(coords::Representation_3D const& cartesians) const override;
    std::tuple<coords::r3, coords::r3, coords::r3, coords::r3>
      der(coords::Representation_3D const& cartesians) const;
    std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override;
    coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
    std::string info(coords::Representation_3D const& cartesians) const override;

	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

    std::size_t index_a_;
    std::size_t index_b_;
    std::size_t index_c_;
    std::size_t index_d_;

    bool operator==(DihedralAngle const&) const;
    
    bool constrained_;
    virtual bool is_constrained() const override {return constrained_;}
  };

  struct OutOfPlane : public DihedralAngle {
    template <typename Atom>
    OutOfPlane(Atom const& outerLeftAtom, Atom const& leftAtom,
      Atom const& rightAtom, Atom const& outerRightAtom)
    : DihedralAngle{ outerLeftAtom, leftAtom, rightAtom, outerRightAtom },
      constrained_{ false /*Config::get().constrained_internals.constrain_out_of_plane_bends*/ }
    {}
    
    using DihedralAngle::DihedralAngle;

    coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
    std::string info(coords::Representation_3D const& cartesians) const override;

	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }
    
    bool constrained_;
    virtual bool is_constrained() const override {return constrained_;}
  };

  struct Translations : public InternalCoordinates::InternalCoordinate {
    virtual ~Translations() = default;

    virtual coords::float_type val(coords::Representation_3D const& cartesians) const override;
    virtual std::string info(coords::Representation_3D const& cartesians) const override;
    virtual std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override;

    std::vector<std::size_t> indices_;

    coords::float_type hessian_guess(coords::Representation_3D const& /*cartesians*/) const override {
      return 0.05;
    }

    bool operator==(Translations const&) const;

	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

    bool constrained_;
    virtual bool is_constrained() const override {return constrained_;}

  protected:
    using CoordinateFunc = scon::_c3::_fc<coords::float_type> (coords::Cartesian_Point::*)() const;

    Translations(std::vector<std::size_t> const& index_vec, CoordinateFunc cf, char letter):
        constrained_{ false /*Config::get().constrained_internals.constrain_translations*/ },
        coord_func_{cf},
        coordinate_letter{letter}
    {
      for (auto index : index_vec) {
        indices_.emplace_back(index - 1u);
      }
    }

  private:
    const CoordinateFunc coord_func_;
    const char coordinate_letter;

    virtual coords::Cartesian_Point size_reciprocal(std::size_t s) const = 0;
  };

  struct TranslationX : Translations {
    explicit TranslationX(std::vector<std::size_t> const& index_vec):
        Translations(index_vec, &coords::Cartesian_Point::x, 'X')
    {}

  private:
    virtual coords::Cartesian_Point size_reciprocal(std::size_t s) const override{
      return coords::Cartesian_Point(1. / static_cast<coords::float_type>(s), 0., 0.);
    }
  };

  struct TranslationY : Translations {
    explicit TranslationY(std::vector<std::size_t> const& index_vec):
        Translations(index_vec, &coords::Cartesian_Point::y, 'Y')
    {}

  private:
    virtual coords::Cartesian_Point size_reciprocal(std::size_t s) const override{
      return coords::Cartesian_Point(0., 1. / static_cast<coords::float_type>(s), 0.);
    }
  };

  struct TranslationZ : Translations {
    explicit TranslationZ(std::vector<std::size_t> const& index_vec):
        Translations(index_vec, &coords::Cartesian_Point::z, 'Z')
    {}

  private:
    virtual coords::Cartesian_Point size_reciprocal(std::size_t s) const override{
      return coords::Cartesian_Point(0., 0., 1. / static_cast<coords::float_type>(s));
    }
  };
  
  template<typename T>
  template<typename T_>
  inline CartesiansForInternalCoordinatesImpl<T>::CartesiansForInternalCoordinatesImpl(T_ && cartesians) 
    : coords::Representation_3D(std::forward<T_>(cartesians)) {}

  template<typename T>
  template<typename T_>
  inline CartesiansForInternalCoordinatesImpl<T>&
    CartesiansForInternalCoordinatesImpl<T>::setCartesianCoordnates(T_ && newCartesianCoordinates) {
    notify();
    return *this = std::forward<T_>(newCartesianCoordinates);
  }

  class Rotator;
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
	//forward declare mathmatrix
    scon::mathmatrix<coords::float_type> const& rot_der_mat(coords::Representation_3D const&);

    Rotations makeRotations();

    bool operator==(Rotator const& other) const;
    
  private:
	Rotator(coords::Representation_3D const& reference, std::vector<std::size_t> const& index_vec);
   
    std::vector<scon::mathmatrix<coords::float_type>> rot_der(coords::Representation_3D const&) const;
    coords::float_type radiusOfGyration(const coords::Representation_3D&);


    std::array<coords::float_type, 3u> storedValuesForRotations;
    std::unique_ptr<scon::mathmatrix<coords::float_type>> storedDerivativesForRotations;
    bool updateStoredValues;
    bool updateStoredDerivatives;

    void registerCartesians(InternalCoordinates::CartesiansForInternalCoordinates & cartesianCoordinates);

    coords::Representation_3D const reference_;
    std::vector<std::size_t> indices_;
    coords::float_type rad_gyr_;
  };

  struct Rotation : public InternalCoordinate {

    virtual coords::float_type val(coords::Representation_3D const& cartesians) const override {
      auto const& returnValues = rotator->valueOfInternalCoordinate(cartesians);
      return returnValues.at(index_);
    }
	virtual std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override;

    virtual coords::float_type hessian_guess(coords::Representation_3D const& /*cartesians*/) const override {
      return 0.05;
    }

    virtual std::string info(coords::Representation_3D const & cartesians) const override {
      std::ostringstream oss;
      oss << "Rotation " << info_letter_<< ": " << val(cartesians) << " | Constrained: " << std::boolalpha << is_constrained();
      return oss.str();
    }

	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

    std::shared_ptr<Rotator> rotator;
    
    bool constrained_;
    virtual bool is_constrained() const override {return constrained_;}

  protected:
    Rotation(std::shared_ptr<Rotator> rotator, size_t index, char info_letter):
        rotator{std::move(rotator)},
        constrained_{false},
        index_{index},
        info_letter_{info_letter}
    {}

  private:
    std::size_t index_;
    char info_letter_;
  };

  struct RotationA : public Rotation {
    explicit RotationA(std::shared_ptr<Rotator> rotator):
        Rotation{std::move(rotator), 0, 'A'}
    {}


    bool operator==(RotationA const& other) const {
      return *rotator.get() == *other.rotator.get();
    }
  };

  struct RotationB : public Rotation {
    explicit RotationB(std::shared_ptr<Rotator> rotator):
        Rotation{std::move(rotator), 1, 'B'}
    {}

    bool operator==(RotationB const& other) const {
      return *rotator.get() == *other.rotator.get();
    }
  };


  struct RotationC : public Rotation {
    explicit RotationC(std::shared_ptr<Rotator> rotator):
        Rotation{std::move(rotator), 2, 'C'}
    {}

    bool operator==(RotationC const& other) const {
      return *rotator.get() == *other.rotator.get();
    }
  };

  template<typename T>
  void CartesiansForInternalCoordinatesImpl<T>::registerObserver(std::shared_ptr<RotatorObserver> const observer) {
    observerList.emplace_back(observer);
  }

  template<typename T>
  void CartesiansForInternalCoordinatesImpl<T>::notify(){
    for (auto const& observer : observerList) {
      observer->update();
    }
  }

  
}

#endif
