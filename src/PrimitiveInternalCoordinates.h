/**
CAST 3
PrimitiveInternalCoordinates.h
Purpose: Definition of primitive Internal Coordinate Systems


@author Julian Erdmannsdörfer, Michael Prem
@version 3.0
*/

#ifndef PRIMITIVE_INTERNAL_COORDINATES_H
#define PRIMITIVE_INTERNAL_COORDINATES_H

#include"InternalCoordinateBase.h"
#include"coords.h"
#include "graph.h"

namespace internals {
  class AppropriateStepFinder;
  class ICDecoratorBase;
  class InternalToCartesianConverter;

  class ICDecoratorBase;

  class PrimitiveInternalCoordinates {
  public:
    PrimitiveInternalCoordinates();
    explicit PrimitiveInternalCoordinates(ICDecoratorBase & decorator);

	virtual ~PrimitiveInternalCoordinates();

    void appendPrimitives(InternalVec && primitives);
    void appendRotators(std::vector<std::shared_ptr<InternalCoordinates::Rotator>> const& rotators);

    virtual std::unique_ptr<AppropriateStepFinder> constructStepFinder(
      InternalToCartesianConverter const& converter,
      scon::mathmatrix<coords::float_type> const& gradients,
      scon::mathmatrix<coords::float_type> const& hessian,
      CartesianType const& /*cartesians*/
    );

    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> primitive_internals;
    //std::vector<std::shared_ptr<InternalCoordinates::Rotator>> rotation_vec_;

    void requestNewBAndG(){
     new_B_matrix = true;
     new_G_matrix = true;
    }

  protected:
    std::vector<std::shared_ptr<InternalCoordinates::Rotator>> registeredRotators;

	std::unique_ptr<scon::mathmatrix<coords::float_type>> B_matrix;
	std::unique_ptr<scon::mathmatrix<coords::float_type>> G_matrix;
	std::unique_ptr<scon::mathmatrix<coords::float_type>> hessian;

    std::vector<std::vector<coords::float_type>> deriv_vec(CartesianType const& cartesians);

    bool new_B_matrix = true;
    bool new_G_matrix = true;
    
  public:
    virtual scon::mathmatrix<coords::float_type> calc(coords::Representation_3D const& xyz) const;//F
    virtual scon::mathmatrix<coords::float_type> calc_diff(coords::Representation_3D const& lhs, coords::Representation_3D const& rhs) const;//F

    virtual scon::mathmatrix<coords::float_type> guess_hessian(CartesianType const&) const;//F
    virtual scon::mathmatrix<coords::float_type>& Bmat(CartesianType const& cartesians);//F
    virtual scon::mathmatrix<coords::float_type>& Gmat(CartesianType const& cartesians);//F
    virtual scon::mathmatrix<coords::float_type> transposeOfBmat(CartesianType const& cartesian);
    virtual scon::mathmatrix<coords::float_type> pseudoInverseOfGmat(CartesianType const& cartesian);
    virtual scon::mathmatrix<coords::float_type> projectorMatrix(CartesianType const& cartesian);

  };

  class InternalToCartesianConverter {
  public:
    InternalToCartesianConverter(PrimitiveInternalCoordinates & internals,
      InternalCoordinates::CartesiansForInternalCoordinates & cartesians) : internalCoordinates{ internals }, cartesianCoordinates{ cartesians } {}
    virtual ~InternalToCartesianConverter();

    scon::mathmatrix<coords::float_type> calculateInternalGradients(scon::mathmatrix<coords::float_type> const&);//Test?

    virtual coords::Representation_3D applyInternalChange(scon::mathmatrix<coords::float_type>) const;
    template<typename XYZ>
    coords::Representation_3D& set_xyz(XYZ&& new_xyz);
    virtual InternalCoordinates::CartesiansForInternalCoordinates const& getCartesianCoordinates() const { return cartesianCoordinates; }
    virtual InternalCoordinates::CartesiansForInternalCoordinates & getCartesianCoordinates() { return cartesianCoordinates; }
    
    std::pair<coords::float_type, coords::float_type> cartesianNormOfOtherStructureAndCurrent(coords::Representation_3D const& otherCartesians) const;//Test
	scon::mathmatrix<coords::float_type> calculateInternalValues()const;
    virtual void reset();
  protected:
    PrimitiveInternalCoordinates & internalCoordinates;
    InternalCoordinates::CartesiansForInternalCoordinates & cartesianCoordinates;

  private:
    coords::Representation_3D& takeCartesianStep(scon::mathmatrix<coords::float_type>&& d_cart);
    void takeCartesianStep(scon::mathmatrix <coords::float_type> && cartesianChange, InternalCoordinates::temporaryCartesian & cartesians) const;
  };

  class RandomNumberForHessianAlteration{
  public:
    RandomNumberForHessianAlteration() : engine{std::random_device()()}, distribution(0.0,1.0){}
    coords::float_type getRandomNumberBetweenZeroAndOne() {
      return distribution(engine);
    }
  private:
    std::mt19937 engine;
    std::uniform_real_distribution<coords::float_type> distribution;
  };

  class AppropriateStepFinder;

  class StepRestrictor {
  public:
	StepRestrictor(std::shared_ptr<scon::mathmatrix<coords::float_type>> step, std::shared_ptr<coords::Representation_3D> cartesians, coords::float_type const target);
	StepRestrictor(StepRestrictor const&) = delete;
	StepRestrictor(StepRestrictor && other);
	virtual ~StepRestrictor();
    virtual coords::float_type operator()(AppropriateStepFinder & finder);

    void registerBestGuess();

    coords::float_type getRestrictedSol() const { return restrictedSol; }
    virtual scon::mathmatrix<coords::float_type> const& getRestrictedStep() const { return *restrictedStep; }
    virtual scon::mathmatrix<coords::float_type> & getRestrictedStep() { return *restrictedStep; }

    void setCartesians(coords::Representation_3D && cartesians) { *correspondingCartesians = std::move(cartesians); }
    coords::Representation_3D const& getCartesians() const { return *correspondingCartesians; }

    void setInitialV0(coords::float_type const initialV0) { v0 = initialV0; }
    
    bool targetIsZero() const { return target == 0.0; }
    
    coords::float_type getTarget() const { return target; }

	void swap(StepRestrictor && other) {
		stepCallbackReference.swap(stepCallbackReference);
		cartesianCallbackReference.swap(cartesianCallbackReference);
		std::swap(target, other.target);
		restrictedStep.swap(other.restrictedStep);
		correspondingCartesians.swap(correspondingCartesians);
		std::swap(restrictedSol, other.restrictedSol);
		std::swap(v0, other.v0);
	}
  
  protected:
    scon::mathmatrix<coords::float_type> alterHessian(scon::mathmatrix<coords::float_type> const & hessian, coords::float_type const alteration) const;
    void randomizeAlteration(std::size_t const step);
	coords::float_type getStepNorm() const;

    std::shared_ptr<scon::mathmatrix<coords::float_type>> stepCallbackReference;
    std::shared_ptr<coords::Representation_3D> cartesianCallbackReference;
    coords::float_type target;
    std::unique_ptr<scon::mathmatrix<coords::float_type>> restrictedStep;
    std::unique_ptr<coords::Representation_3D> correspondingCartesians;
    coords::float_type restrictedSol, v0;
  };

  class StepRestrictorFactory {
  public:
    StepRestrictorFactory(AppropriateStepFinder & finder);
	virtual ~StepRestrictorFactory();
    StepRestrictor makeStepRestrictor(coords::float_type const target);

  private:
    std::shared_ptr<scon::mathmatrix<coords::float_type>> finalStep;
    std::shared_ptr<coords::Representation_3D> finalCartesians;
  };

  class InternalToCartesianStep {
  public:
    InternalToCartesianStep(AppropriateStepFinder & finder, coords::float_type const trustRadius)
      : finder{ finder }, trustRadius {trustRadius} {}
    virtual ~InternalToCartesianStep();

    virtual coords::float_type operator()(StepRestrictor & restrictor);
    
  protected:
    AppropriateStepFinder & finder;
    coords::float_type trustRadius;
  };

  class BrentsMethod {
  public:
    BrentsMethod(AppropriateStepFinder & finder, coords::float_type const leftLimit, coords::float_type const rightLimit, coords::float_type const trustStep)
      : finder{ finder }, leftLimit{ leftLimit }, middle{ 0.0 }, oldMiddle{ 0.0 }, rightLimit{ rightLimit }, result{ 0.0 },
      trustStep{ trustStep }, threshold{ 0.1 }, delta{ 1.e-6 }, bisectionWasUsed{ true } {}

    BrentsMethod(AppropriateStepFinder & finder, coords::float_type const leftLimit, coords::float_type const rightLimit, coords::float_type const trustStep, coords::float_type const cartesianNorm)
      : finder{ finder }, leftLimit{ leftLimit }, middle{ 0.0 }, oldMiddle{ 0.0 }, rightLimit{ rightLimit }, result{ 0.0 },
      valueLeft{ -trustStep }, valueRight{ cartesianNorm - trustStep}, trustStep{ trustStep }, threshold{ 0.1 }, delta{ 1.e-6 }, bisectionWasUsed{ true } {}

	virtual ~BrentsMethod();

    coords::float_type operator()(InternalToCartesianStep & internalToCartesianStep);
  protected:
    bool useBisection()const;

    AppropriateStepFinder & finder;
    coords::float_type leftLimit, middle, oldMiddle, rightLimit, result, valueLeft, valueRight;
    coords::float_type const trustStep;
    coords::float_type const threshold;
    coords::float_type const delta;
    bool bisectionWasUsed;
  };
  

  struct GradientsAndHessians;

  class AppropriateStepFinder {
  public:
	AppropriateStepFinder(InternalToCartesianConverter const& converter, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian);


	std::unique_ptr<GradientsAndHessians> matrices;


    virtual void appropriateStep(coords::float_type const trustRadius);

    virtual coords::float_type getDeltaYPrime(scon::mathmatrix<coords::float_type> const& internalStep) const;
    virtual coords::float_type getSol(scon::mathmatrix<coords::float_type> const& internalStep) const;
    virtual coords::float_type getSolBestStep() const;//Test

    virtual scon::mathmatrix<coords::float_type> getInternalStep() const;
    virtual scon::mathmatrix<coords::float_type> getInternalStep(scon::mathmatrix<coords::float_type> const& hessian) const;

    coords::Representation_3D & getCartesians() { return *bestCartesiansSoFar; }

    StepRestrictor generateStepRestrictor(coords::float_type const target);

    virtual coords::float_type applyInternalChangeAndGetNorm(StepRestrictor & internalStep);
    virtual coords::float_type applyInternalChangeAndGetNorm(scon::mathmatrix<coords::float_type> const& internalStep);

    virtual scon::mathmatrix<coords::float_type> alterHessian(coords::float_type const alteration) const;


    scon::mathmatrix<coords::float_type> && extractBestStep() { return std::move(*bestStepSoFar); }
    coords::Representation_3D && extractCartesians() { return std::move(*bestCartesiansSoFar); }
    scon::mathmatrix<coords::float_type> const& getBestStep() const {return *bestStepSoFar;}

	virtual ~AppropriateStepFinder();

  protected:
    //Constructor for Testclass
	AppropriateStepFinder(InternalToCartesianConverter const& converter, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian, scon::mathmatrix<coords::float_type> && invertedHessian);

    friend class StepRestrictorFactory;

	std::shared_ptr<scon::mathmatrix<coords::float_type>> getAddressOfInternalStep() { return bestStepSoFar; }
	std::shared_ptr<coords::Representation_3D> getAddressOfCartesians() { return bestCartesiansSoFar; }

    InternalToCartesianConverter const& converter;
    std::shared_ptr<scon::mathmatrix<coords::float_type>> bestStepSoFar;
	std::shared_ptr<coords::Representation_3D> bestCartesiansSoFar;
    StepRestrictorFactory stepRestrictorFactory;
  };


  template<typename XYZ>
  coords::Representation_3D& InternalToCartesianConverter::set_xyz(XYZ&& new_xyz) {
    internalCoordinates.requestNewBAndG();
    return cartesianCoordinates.setCartesianCoordnates(std::forward<XYZ>(new_xyz));
  }
}

#endif
