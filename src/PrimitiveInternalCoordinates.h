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
#include "ic_core.h"
#include "graph.h"

namespace internals {
  class PrimitiveInternalCoordinates : public InternalCoordinatesBase, public std::enable_shared_from_this<PrimitiveInternalCoordinates> {
  public:
    PrimitiveInternalCoordinates() = default;
    virtual ~PrimitiveInternalCoordinates() = default;
    
    
    virtual void buildCoordinates(CartesianType & /*cartesians*/, BondGraph const& /*graph*/, IndexVec const& /*indexVec*/) override{} /// We are not building any coordinates here. Everything is done by the decorators.
    void appendCoordinates(std::shared_ptr<InternalCoordinateAppenderInterface> appender) override;
    void appendPrimitives(InternalVec && primitives);
    void appendRotators(std::vector<std::shared_ptr<InternalCoordinates::Rotator>> const& rotators);

    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> primitive_internals;
    //std::vector<std::shared_ptr<InternalCoordinates::Rotator>> rotation_vec_;

    void requestNewBAndG(){
     new_B_matrix = true;
     new_G_matrix = true;
    }

  protected:
    std::vector<std::shared_ptr<InternalCoordinates::Rotator>> registeredRotators;

    scon::mathmatrix<coords::float_type> B_matrix;
    scon::mathmatrix<coords::float_type> G_matrix;
    scon::mathmatrix<coords::float_type> hessian;

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

  };

  class InternalToCartesianConverter {
  public:
    InternalToCartesianConverter(PrimitiveInternalCoordinates & internals,
      InternalCoordinates::CartesiansForInternalCoordinates & cartesians) : internalCoordinates{ internals }, cartesianCoordinates{ cartesians } {}
    virtual ~InternalToCartesianConverter() = default;

    scon::mathmatrix<coords::float_type> calculateInternalGradients(scon::mathmatrix<coords::float_type> const&);//Test?

    virtual coords::Representation_3D applyInternalChange(scon::mathmatrix<coords::float_type>) const;
    template<typename XYZ>
    coords::Representation_3D& set_xyz(XYZ&& new_xyz);
    virtual InternalCoordinates::CartesiansForInternalCoordinates const& getCartesianCoordinates() const { return cartesianCoordinates; }
    virtual InternalCoordinates::CartesiansForInternalCoordinates & getCartesianCoordinates() { return cartesianCoordinates; }
    
    std::pair<coords::float_type, coords::float_type> cartesianNormOfOtherStructureAndCurrent(coords::Representation_3D const& otherCartesians) const;//Test
    scon::mathmatrix<coords::float_type> calculateInternalValues()const{
      return internalCoordinates.calc(cartesianCoordinates);
    }
    virtual void reset();
  protected:
    PrimitiveInternalCoordinates & internalCoordinates;
    InternalCoordinates::CartesiansForInternalCoordinates & cartesianCoordinates;

  private:
    template<typename Dcart>
    coords::Representation_3D& takeCartesianStep(Dcart&& d_cart);
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
    StepRestrictor(scon::mathmatrix<coords::float_type> * step, coords::Representation_3D * cartesians, coords::float_type const target) : stepCallbackReference{ step }, cartesianCallbackReference{ cartesians }, target{ target }, restrictedStep{}, correspondingCartesians{}, restrictedSol { 0.0 }, v0{ 0.0 } {}
    virtual ~StepRestrictor() = default;
    virtual coords::float_type operator()(AppropriateStepFinder & finder);

    void registerBestGuess();

    coords::float_type getRestrictedSol() const { return restrictedSol; }
    virtual scon::mathmatrix<coords::float_type> const& getRestrictedStep() const { return restrictedStep; }
    virtual scon::mathmatrix<coords::float_type> & getRestrictedStep() { return restrictedStep; }

    void setCartesians(coords::Representation_3D && cartesians) { correspondingCartesians = std::move(cartesians); }
    coords::Representation_3D const& getCartesians() const { return correspondingCartesians; }

    void setInitialV0(coords::float_type const initialV0) { v0 = initialV0; }
    
    bool targetIsZero() const { return target == 0.0; }
    
    coords::float_type getTarget() const { return target; }
  
  protected:
    scon::mathmatrix<coords::float_type> alterHessian(scon::mathmatrix<coords::float_type> const & hessian, coords::float_type const alteration) const;
    void randomizeAlteration(std::size_t const step);
    coords::float_type getStepNorm() const { return restrictedStep.norm(); }

    scon::mathmatrix<coords::float_type> * stepCallbackReference;//TODO make these pointers to shared pointer
    coords::Representation_3D * cartesianCallbackReference;//TODO make these pointers to shared pointer
    coords::float_type target;
    scon::mathmatrix<coords::float_type> restrictedStep;
    coords::Representation_3D correspondingCartesians;
    coords::float_type restrictedSol, v0;
  };

  class StepRestrictorFactory {
  public:
    StepRestrictorFactory(AppropriateStepFinder & finder);

    StepRestrictor makeStepRestrictor(coords::float_type const target);

  private:
    scon::mathmatrix<coords::float_type> * const finalStep;
    coords::Representation_3D * const finalCartesians;
  };

  class InternalToCartesianStep {
  public:
    InternalToCartesianStep(AppropriateStepFinder & finder, coords::float_type const trustRadius)
      : finder{ finder }, trustRadius {trustRadius} {}
    virtual ~InternalToCartesianStep() = default;

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
      trustStep{ trustStep }, threshold{ 0.1 }, delta{ 1.e-6 }, bisectionWasUsed{ true }, valueLeft{ -trustStep }, valueRight{ cartesianNorm - trustStep} {}

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
  
  class AppropriateStepFinder {
  public:
    AppropriateStepFinder(InternalToCartesianConverter const& converter, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian) : 
      converter{ converter }, gradients { gradients }, hessian{ hessian }, inverseHessian{ hessian.pinv() }, bestStepSoFar{}, stepRestrictorFactory{ *this } {}

    scon::mathmatrix<coords::float_type> const& gradients;
    scon::mathmatrix<coords::float_type> const& hessian;
    scon::mathmatrix<coords::float_type> inverseHessian;

    virtual void appropriateStep(coords::float_type const trustRadius);

    virtual coords::float_type getDeltaYPrime(scon::mathmatrix<coords::float_type> const& internalStep) const;
    virtual coords::float_type getSol(scon::mathmatrix<coords::float_type> const& internalStep) const;
    virtual coords::float_type getSolBestStep() const;//Test

    virtual scon::mathmatrix<coords::float_type> getInternalStep() const;
    virtual scon::mathmatrix<coords::float_type> getInternalStep(scon::mathmatrix<coords::float_type> const& hessian) const;

    coords::Representation_3D & getCartesians() { return bestCartesiansSoFar; }

    StepRestrictor generateStepRestrictor(coords::float_type const target);

    virtual coords::float_type applyInternalChangeAndGetNorm(StepRestrictor & internalStep);
    virtual coords::float_type applyInternalChangeAndGetNorm(scon::mathmatrix<coords::float_type> const& internalStep);

    virtual scon::mathmatrix<coords::float_type> alterHessian(coords::float_type const alteration) const;


    scon::mathmatrix<coords::float_type> && extractBestStep() { return std::move(bestStepSoFar); }
    coords::Representation_3D && extractCartesians() { return std::move(bestCartesiansSoFar); }
    scon::mathmatrix<coords::float_type> const& getBestStep() const {return bestStepSoFar;}

  protected:
    //Constructor for Testclass
    AppropriateStepFinder(InternalToCartesianConverter const& converter, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian, scon::mathmatrix<coords::float_type> && invertedHessian) :
      converter{ converter }, gradients{ gradients }, hessian{ hessian }, inverseHessian{ std::move(invertedHessian) }, bestStepSoFar{}, stepRestrictorFactory{ *this } {}

    friend class StepRestrictorFactory;

    scon::mathmatrix<coords::float_type> * getAddressOfInternalStep() { return std::addressof(bestStepSoFar); }
    coords::Representation_3D * getAddressOfCartesians() { return std::addressof(bestCartesiansSoFar); }

    InternalToCartesianConverter const& converter;
    scon::mathmatrix<coords::float_type> bestStepSoFar;
    coords::Representation_3D bestCartesiansSoFar;
    StepRestrictorFactory stepRestrictorFactory;
  };

  template<typename Dcart>
  coords::Representation_3D& InternalToCartesianConverter::takeCartesianStep(Dcart&& d_cart) {
    auto d_cart_rep3D = ic_util::matToRep3D(std::forward<Dcart>(d_cart));
    return set_xyz(cartesianCoordinates + d_cart_rep3D);
  }

  template<typename XYZ>
  coords::Representation_3D& InternalToCartesianConverter::set_xyz(XYZ&& new_xyz) {
    internalCoordinates.requestNewBAndG();
    return cartesianCoordinates.setCartesianCoordnates(std::forward<XYZ>(new_xyz));
  }
}

#endif
