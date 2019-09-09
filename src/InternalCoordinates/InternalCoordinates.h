#ifndef H_INTERNAL_COORDINATES
#define H_INTERNAL_COORDINATES

#include<array>

#include"InternalCoordinatesAliases.h"
#include "BondGraph/ElementInformations.h"
#include "../coords.h"

class AbstractConstraintManager;

namespace ic_util {
	enum class period;
}

namespace InternalCoordinates {

  struct InternalCoordinate {
	virtual coords::float_type val(coords::Representation_3D const& cartesians) const = 0;
	virtual coords::float_type difference(coords::Representation_3D const& newCoordinates, coords::Representation_3D const& oldCoordinates) const = 0;
	virtual std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const = 0;
	virtual coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const = 0;
	virtual std::string info(coords::Representation_3D const & cartesians) const = 0;
	virtual bool hasIndices(std::vector<std::size_t> const& indices) const = 0;
	virtual std::vector<std::size_t> getIndices() const = 0;
	virtual void makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) = 0;
	virtual void makeConstrained() = 0;
	virtual void releaseConstraint() = 0;
	virtual bool is_constrained() const = 0;
	virtual ~InternalCoordinate() = 0;
  };

	InternalCoordinate::~InternalCoordinate() = default;

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

  class CartesiansForInternalCoordinates {
  public:


		CartesiansForInternalCoordinates(CartesiansForInternalCoordinates && cartesians);
	  CartesiansForInternalCoordinates(CartesiansForInternalCoordinates const& cartesians);

	  CartesiansForInternalCoordinates(coords::Representation_3D && cartesians);
	  CartesiansForInternalCoordinates(coords::Representation_3D const& cartesians);

	  CartesiansForInternalCoordinates & operator=(CartesiansForInternalCoordinates const& cartesians);
	  CartesiansForInternalCoordinates & operator=(CartesiansForInternalCoordinates && cartesians);

		CartesiansForInternalCoordinates() = default;
		virtual ~CartesiansForInternalCoordinates() = default;

		friend coords::Representation_3D operator+(CartesiansForInternalCoordinates const& lhs, coords::Representation_3D const& rhs);

		coords::Cartesian_Point const& at(std::size_t const i) const;
		coords::Cartesian_Point & at(std::size_t const i);

		coords::float_type getInternalValue(InternalCoordinate const& in) const;
		std::vector<coords::float_type> getInternalDerivativeVector(InternalCoordinate const& in) const;
		coords::float_type getInternalHessianGuess(InternalCoordinate const& in) const;
		coords::float_type getInternalDifference(CartesiansForInternalCoordinates const& other, InternalCoordinate const& in) const;

		std::pair<coords::float_type, coords::float_type> displacementRmsValAndMaxTwoStructures(coords::Representation_3D const& other) const;

		std::pair<coords::float_type, coords::float_type> displacementRmsValAndMaxTwoStructures(CartesiansForInternalCoordinates const& other) const;

		coords::Representation_3D toAngstrom() const;

		void registerObserver(std::shared_ptr<RotatorObserver> const observer);

		void setCartesianCoordnates(coords::Representation_3D const& newCartesianCoordinates);
		void setCartesianCoordnates(coords::Representation_3D&& newCartesianCoordinates);

		void reset();

  protected:
    friend class temporaryCartesian;
    std::vector<std::shared_ptr<AbstractGeometryObserver>> observerList;
	coords::Representation_3D coordinates;
  private:
	template<typename T_>
	void setCartesianCoordnatesIntern(T_ && newCartesianCoordinates) {
	  notify();
	  coordinates = std::forward<T_>(newCartesianCoordinates);
	}
    void notify();
  };

  class temporaryCartesian{
  public:
	  temporaryCartesian(CartesiansForInternalCoordinates & cartesians) : coordinates{ cartesians.coordinates }, coordinatesPtr{ &cartesians } {}
	coords::Container<coords::Cartesian_Point> coordinates;
	void stolenNotify() {
		coordinatesPtr->notify();
	}
  private:
	  CartesiansForInternalCoordinates * coordinatesPtr;
  };

  inline coords::Representation_3D sliceCartesianCoordinates(CartesiansForInternalCoordinates const& cartesians, std::vector<std::size_t> const& indexVector) {
    coords::Representation_3D slicedCoordinates;
    for (auto const& index : indexVector) {
      slicedCoordinates.emplace_back(cartesians.at(index - 1));
    }
    return slicedCoordinates;
  }

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
	coords::float_type difference(coords::Representation_3D const& newCoordinates, coords::Representation_3D const& oldCoordinates) const override;
    std::pair<coords::r3, coords::r3> der(coords::Representation_3D const& cartesians) const;
    std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override;
    coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
    std::string info(coords::Representation_3D const& cartesians) const override;

		virtual bool hasIndices(std::vector<std::size_t> const& indices) const override;
		virtual std::vector<std::size_t> getIndices() const override;
		virtual void makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) override;
	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

    bool operator==(BondDistance const&) const;
    
    bool constrained_;
    virtual bool is_constrained() const override {return constrained_;}

  private:
    bool bothElementsInPeriodOne(period const atomA, period const atomB)const;
    bool oneElementInPeriodOneTheOtherInPeriodTwo(period const atomA, period const atomB)const;
    bool oneElementInPeriodOneTheOtherInPeriodThree(period const atomA, period const atomB)const;
    bool bothElementsInPeriodTwo(period const atomA, period const atomB)const;
    bool oneElementInPeriodTwoTheOtherInPeriodThree(period const atomA, period const atomB)const;
  };

  struct BondAngle : InternalCoordinate {
    template<typename Atom>
    BondAngle(Atom const& leftAtom, Atom const& middleAtom, Atom const& rightAtom)
      : index_a_{ leftAtom.atom_serial - 1u }, index_b_{ middleAtom.atom_serial - 1u },
      index_c_{ rightAtom.atom_serial - 1u }, elem_a_{ leftAtom.element }, elem_b_{ middleAtom.element },
      elem_c_{ rightAtom.element },
      constrained_{ false }
    {}

    std::size_t index_a_;
    std::size_t index_b_;
    std::size_t index_c_;
    std::string elem_a_;
    std::string elem_b_;
    std::string elem_c_;

    coords::float_type val(coords::Representation_3D const& cartesians) const override;
	coords::float_type difference(coords::Representation_3D const& newCoordinates, coords::Representation_3D const& oldCoordinates) const override;
    std::tuple<coords::r3, coords::r3, coords::r3> der(coords::Representation_3D const& cartesians) const;
    std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override;
    coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
    std::string info(coords::Representation_3D const& cartesians) const override;

		virtual bool hasIndices(std::vector<std::size_t> const& indices) const override;
		virtual std::vector<std::size_t> getIndices() const override;
		virtual void makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) override;
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
      constrained_{ false }{}
    virtual ~DihedralAngle() = default;

    coords::float_type val(coords::Representation_3D const& cartesians) const override;
	coords::float_type difference(coords::Representation_3D const& newCoordinates, coords::Representation_3D const& oldCoordinates) const override;
    std::tuple<coords::r3, coords::r3, coords::r3, coords::r3>
      der(coords::Representation_3D const& cartesians) const;
    std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override;
    coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
    std::string info(coords::Representation_3D const& cartesians) const override;

		virtual bool hasIndices(std::vector<std::size_t> const& indices) const override;
		virtual std::vector<std::size_t> getIndices() const override;
		virtual void makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) override;
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
    : DihedralAngle{ outerLeftAtom, leftAtom, rightAtom, outerRightAtom }
    {}
    
    using DihedralAngle::DihedralAngle;

    coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
    std::string info(coords::Representation_3D const& cartesians) const override;

    virtual bool is_constrained() const override {return constrained_;}
  };

  struct Translations : public InternalCoordinates::InternalCoordinate {

		enum class Direction: int {X, Y, Z};


    virtual ~Translations() = default;

    virtual coords::float_type val(coords::Representation_3D const& cartesians) const override;
	coords::float_type difference(coords::Representation_3D const& newCoordinates, coords::Representation_3D const& oldCoordinates) const override;
    virtual std::string info(coords::Representation_3D const& cartesians) const override;
    virtual std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override;

    std::vector<std::size_t> indices_;

    coords::float_type hessian_guess(coords::Representation_3D const& /*cartesians*/) const override {
      return 0.05;
    }

    bool operator==(Translations const&) const;

		virtual bool hasIndices(std::vector<std::size_t> const& indices) const override;
		virtual std::vector<std::size_t> getIndices() const override;
		virtual void makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) override;
	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

    bool constrained_;
    virtual bool is_constrained() const override {return constrained_;}

		virtual Direction getDirection() const = 0;

  protected:
    Translations(std::vector<std::size_t> const& index_vec):
        constrained_{ false }
    {
      for (auto index : index_vec) {
        indices_.emplace_back(index - 1u);
      }
    }
	virtual coords::float_type coord_func(coords::Cartesian_Point const& cp) const = 0;
	virtual char coordinate_letter() const = 0;
	virtual coords::Cartesian_Point size_reciprocal() const = 0;
  };

  struct TranslationX : Translations {
    explicit TranslationX(std::vector<std::size_t> const& index_vec):
        Translations(index_vec)
    {}

  protected:
	coords::float_type coord_func(coords::Cartesian_Point const& cp) const override {
      return cp.x();
	}

	Direction getDirection() const override {return Direction::X;}

	char coordinate_letter() const override { return 'X'; }

	coords::Cartesian_Point size_reciprocal() const override { return coords::Cartesian_Point{ 1./coords::float_type(indices_.size()), 0., 0.}; }
  };

  struct TranslationY : Translations {
    explicit TranslationY(std::vector<std::size_t> const& index_vec):
        Translations(index_vec)
    {}

    protected:
	  coords::float_type coord_func(coords::Cartesian_Point const& cp) const override {
		  return cp.y();
	  }

		Direction getDirection() const override { return Direction::Y; }

	  char coordinate_letter() const override { return 'Y'; }

	  coords::Cartesian_Point size_reciprocal() const override { return coords::Cartesian_Point{ 0., 1. / coords::float_type(indices_.size()), 0. }; }
  };

  struct TranslationZ : Translations {
	explicit TranslationZ(std::vector<std::size_t> const& index_vec) :
		Translations(index_vec)
	{}
  protected:
	coords::float_type coord_func(coords::Cartesian_Point const& cp) const override {
		return cp.z();
	}

	Direction getDirection() const override { return Direction::Z; }

	char coordinate_letter() const override { return 'Z'; }

	coords::Cartesian_Point size_reciprocal() const override { return coords::Cartesian_Point{ 0., 0., 1. / coords::float_type(indices_.size()) }; }
  };
  
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
    scon::mathmatrix<coords::float_type> const& rot_der_mat(coords::Representation_3D const&);

    Rotations makeRotations();
	virtual ~Rotator();

    bool operator==(Rotator const& other) const;
    
  private:

	friend struct Rotation;

	std::array<coords::float_type, 3u> calculateValueOfInternalCoordinate(coords::Representation_3D const& newXyz) const;
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

		enum class Direction : int { A,B,C };

    virtual coords::float_type val(coords::Representation_3D const& cartesians) const override {
      auto const& returnValues = rotator->valueOfInternalCoordinate(cartesians);
      return returnValues.at(index());
    }
	coords::float_type difference(coords::Representation_3D const& newCoordinates, coords::Representation_3D const& oldCoordinates) const override {
		auto const previousVals = rotator->calculateValueOfInternalCoordinate(oldCoordinates);
		return val(newCoordinates) - previousVals.at(index());
	}
	virtual std::vector<coords::float_type> der_vec(coords::Representation_3D const& cartesians) const override;

    virtual coords::float_type hessian_guess(coords::Representation_3D const& /*cartesians*/) const override {
      return 0.05;
    }

    virtual std::string info(coords::Representation_3D const & cartesians) const override {
      std::ostringstream oss;
      oss << "Rotation " << name() << ": " << val(cartesians) << " | Constrained: " << std::boolalpha << is_constrained();
      return oss.str();
    }

		virtual bool hasIndices(std::vector<std::size_t> const& indices) const override;
		virtual std::vector<std::size_t> getIndices() const override;
		virtual void makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) override;
	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

    std::shared_ptr<Rotator> rotator;
    
    bool constrained_;
    virtual bool is_constrained() const override {return constrained_;}

		virtual Direction getDirection() const = 0;

  protected:
    Rotation(std::shared_ptr<Rotator> rotator):
        rotator{std::move(rotator)},
        constrained_{false}
    {}

	virtual std::size_t index() const = 0;
	virtual char name() const = 0;
  };

  struct RotationA : public Rotation {
    explicit RotationA(std::shared_ptr<Rotator> rotator):
        Rotation{std::move(rotator)}
    {}


    bool operator==(RotationA const& other) const {
      return *rotator.get() == *other.rotator.get();
    }
		Direction getDirection() const override {return Direction::A; }
  protected:
	std::size_t index() const override { return 0u; }
	char name() const override { return 'A'; }
  };

  struct RotationB : public Rotation {
    explicit RotationB(std::shared_ptr<Rotator> rotator):
        Rotation{std::move(rotator)}
    {}

    bool operator==(RotationB const& other) const {
      return *rotator.get() == *other.rotator.get();
    }

		Direction getDirection() const override { return Direction::B; }
  protected:
	  std::size_t index() const override { return 1u; }
	  char name() const override { return 'B'; }
  };


  struct RotationC : public Rotation {
    explicit RotationC(std::shared_ptr<Rotator> rotator):
        Rotation{std::move(rotator)}
    {}

    bool operator==(RotationC const& other) const {
      return *rotator.get() == *other.rotator.get();
    }

		Direction getDirection() const override { return Direction::C; }
  protected:
	  std::size_t index() const override { return 2u; }
	  char name() const override { return 'C'; }
  };


	class AbstractInternalCoordinatesBuilder {
	public:
		using TwoInternals = std::pair<std::unique_ptr<InternalCoordinate>, std::unique_ptr<InternalCoordinate>>;
		using ThreeInternals = std::tuple<std::unique_ptr<InternalCoordinate>, std::unique_ptr<InternalCoordinate>, std::unique_ptr<InternalCoordinate>>;
		virtual std::unique_ptr<InternalCoordinate> buildBondDistance(std::size_t const, std::size_t const) const = 0;
		virtual std::unique_ptr<InternalCoordinate> buildBondAngle(std::size_t const, std::size_t const, std::size_t const) const = 0;
		virtual std::unique_ptr<InternalCoordinate> buildDihedralAngle(std::size_t const, std::size_t const, std::size_t const, std::size_t const) const = 0;
		virtual std::unique_ptr<InternalCoordinate> buildTranslationX(std::vector<std::size_t> const&) const = 0;
		virtual std::unique_ptr<InternalCoordinate> buildTranslationY(std::vector<std::size_t> const&) const = 0;
		virtual std::unique_ptr<InternalCoordinate> buildTranslationZ(std::vector<std::size_t> const&) const = 0;
		virtual TwoInternals buildTranslationXY(std::vector<std::size_t> const&) const = 0;
		virtual TwoInternals buildTranslationXZ(std::vector<std::size_t> const&) const = 0;
		virtual TwoInternals buildTranslationYZ(std::vector<std::size_t> const&) const = 0;
		virtual ThreeInternals buildTranslationXYZ(std::vector<std::size_t> const&) const = 0;
		virtual std::unique_ptr<InternalCoordinate> buildRotationA(std::vector<std::size_t> const&) = 0;
		virtual std::unique_ptr<InternalCoordinate> buildRotationB(std::vector<std::size_t> const&) = 0;
		virtual std::unique_ptr<InternalCoordinate> buildRotationC(std::vector<std::size_t> const&) = 0;
		virtual TwoInternals buildRotationAB(std::vector<std::size_t> const&) = 0;
		virtual TwoInternals buildRotationAC(std::vector<std::size_t> const&) = 0;
		virtual TwoInternals buildRotationBC(std::vector<std::size_t> const&) = 0;
		virtual ThreeInternals buildRotationABC(std::vector<std::size_t> const&) = 0;
		virtual ~AbstractInternalCoordinatesBuilder() = 0;
	};

	AbstractInternalCoordinatesBuilder::~AbstractInternalCoordinatesBuilder() = default;

	class InternalCoordinatesBuilder : public AbstractInternalCoordinatesBuilder {
	public:
		InternalCoordinatesBuilder(ic_util::BondGraph const& graph, InternalCoordinates::CartesiansForInternalCoordinates & coordinates, std::vector<std::shared_ptr<InternalCoordinates::Rotator>> & rotators);
		virtual std::unique_ptr<InternalCoordinate> buildBondDistance(std::size_t const, std::size_t const) const override;
		virtual std::unique_ptr<InternalCoordinate> buildBondAngle(std::size_t const, std::size_t const, std::size_t const) const override;
		virtual std::unique_ptr<InternalCoordinate> buildDihedralAngle(std::size_t const, std::size_t const, std::size_t const, std::size_t const) const override;
		virtual std::unique_ptr<InternalCoordinate> buildTranslationX(std::vector<std::size_t> const&) const override;
		virtual std::unique_ptr<InternalCoordinate> buildTranslationY(std::vector<std::size_t> const&) const override;
		virtual std::unique_ptr<InternalCoordinate> buildTranslationZ(std::vector<std::size_t> const&) const override;
		virtual TwoInternals buildTranslationXY(std::vector<std::size_t> const&) const override;
		virtual TwoInternals buildTranslationXZ(std::vector<std::size_t> const&) const override;
		virtual TwoInternals buildTranslationYZ(std::vector<std::size_t> const&) const override;
		virtual ThreeInternals buildTranslationXYZ(std::vector<std::size_t> const&) const override;
		virtual std::unique_ptr<InternalCoordinate> buildRotationA(std::vector<std::size_t> const&) override;
		virtual std::unique_ptr<InternalCoordinate> buildRotationB(std::vector<std::size_t> const&) override;
		virtual std::unique_ptr<InternalCoordinate> buildRotationC(std::vector<std::size_t> const&) override;
		virtual TwoInternals buildRotationAB(std::vector<std::size_t> const&) override;
		virtual TwoInternals buildRotationAC(std::vector<std::size_t> const&) override;
		virtual TwoInternals buildRotationBC(std::vector<std::size_t> const&) override;
		virtual ThreeInternals buildRotationABC(std::vector<std::size_t> const&) override;
	private:
		ic_util::BondGraph const& bondGraph;
		InternalCoordinates::CartesiansForInternalCoordinates & cartesians;
		std::vector<std::shared_ptr<InternalCoordinates::Rotator>> & rotatorContainer;
		std::size_t numberOfAtoms;
	};
}

#endif
