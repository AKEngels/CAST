#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_ROTATIONS_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_ROTATIONS_H_

#include "InternalCoordinate.h"
#include "RotatorListener.h"

#include<memory>
#include<array>

namespace internals {

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


class Rotator : public RotatorListener, public std::enable_shared_from_this<Rotator> {
public:
	static std::shared_ptr<Rotator> buildRotator(InternalCoordinates::CartesiansForInternalCoordinates & cartesians, std::vector<std::size_t> const& indexVector) {
		//The ctor should keep being private. Thus make_shared connot be used. If someone has a better idea, go ahead and refactor :)
		auto newInstance = std::shared_ptr<Rotator>(new Rotator(sliceCartesianCoordinates(cartesians, indexVector), indexVector));
		newInstance->registerCartesians(cartesians);
		return newInstance;
	}

	void setAllFlag()override { updateStoredValues = updateStoredDerivatives = true; }
	void requestNewValueEvaluation() { updateStoredValues = true; }
	bool areValuesUpToDate() { return !updateStoredValues; }
	bool areDerivativesUpToDate() { return !updateStoredDerivatives; }

	std::array<float_type, 3u> const& valueOfInternalCoordinate(const coords::Representation_3D&);
	scon::mathmatrix<float_type> const& rot_der_mat(coords::Representation_3D const&);

	Rotations makeRotations();
	virtual ~Rotator();

	bool operator==(Rotator const& other) const;

private:

	friend struct Rotation;

	std::array<float_type, 3u> calculateValueOfInternalCoordinate(coords::Representation_3D const& newXyz) const;
	Rotator(coords::Representation_3D const& reference, std::vector<std::size_t> const& index_vec);

	std::vector<scon::mathmatrix<float_type>> rot_der(coords::Representation_3D const&) const;
	float_type radiusOfGyration(const coords::Representation_3D&);


	std::array<float_type, 3u> storedValuesForRotations;
	std::unique_ptr<scon::mathmatrix<float_type>> storedDerivativesForRotations;
	bool updateStoredValues;
	bool updateStoredDerivatives;

	void registerCartesians(InternalCoordinates::CartesiansForInternalCoordinates & cartesianCoordinates);

	std::unique_ptr<coords::Representation_3D> reference_;
	std::vector<std::size_t> indices_;
	float_type rad_gyr_;
};

struct Rotation : public InternalCoordinate {

	enum class Direction : int { A, B, C };

	virtual float_type val(scon::mathmatrix<float_type> const& cartesians) const override;
	virtual float_type difference(scon::mathmatrix<float_type> const& newCoordinates, scon::mathmatrix<float_type> const& oldCoordinates) const override;
	virtual scon::mathmatrix<float_type> der_vec(scon::mathmatrix<float_type> const& cartesians) const override;

	virtual  float_type hessian_guess(scon::mathmatrix<float_type> const& /*cartesians*/) const override {
		return 0.05;
	}

	virtual std::string info(scon::mathmatrix<float_type> const & cartesians) const override;

	bool hasIndices(std::vector<std::size_t> const& indices) const;
	std::vector<std::size_t> getIndices() const;
	virtual void makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) override;
	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

	std::shared_ptr<Rotator> rotator;

	bool constrained_;
	virtual bool is_constrained() const override { return constrained_; }

	virtual Direction getDirection() const = 0;

protected:
	Rotation(std::shared_ptr<Rotator> rotator) :
		rotator{ std::move(rotator) },
		constrained_{ false }
	{}

	virtual std::size_t index() const = 0;
	virtual char name() const = 0;
};

struct RotationA : public Rotation {
	explicit RotationA(std::shared_ptr<Rotator> rotator) :
		Rotation{ std::move(rotator) }
	{}


	bool operator==(RotationA const& other) const {
		return *rotator.get() == *other.rotator.get();
	}
	Direction getDirection() const override { return Direction::A; }
protected:
	std::size_t index() const override { return 0u; }
	char name() const override { return 'A'; }
};

struct RotationB : public Rotation {
	explicit RotationB(std::shared_ptr<Rotator> rotator) :
		Rotation{ std::move(rotator) }
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
	explicit RotationC(std::shared_ptr<Rotator> rotator) :
		Rotation{ std::move(rotator) }
	{}

	bool operator==(RotationC const& other) const {
		return *rotator.get() == *other.rotator.get();
	}

	Direction getDirection() const override { return Direction::C; }
protected:
	std::size_t index() const override { return 2u; }
	char name() const override { return 'C'; }
};

}

#endif // CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_ROTATIONS_H_