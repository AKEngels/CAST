#ifdef GOOGLE_MOCK

#ifndef PRIMITIVE_INTERNALS_TEST_H
#define PRIMITIVE_INTERNALS_TEST_H

#include<gtest/gtest.h>
#include<gmock/gmock.h>
#include<vector>
#include<utility>

#include"../../graph.h"
#include"../../coords.h"
#include"../../PrimitiveInternalCoordinates.h"
#include"../../TranslationRotationInternalCoordinates.h"
#include"../../InternalCoordinates.h"
#include"ExpectedValuesForInternalCoordinatesTest.h"

class MockPrimitiveInternals : public internals::PrimitiveInternalCoordinates {
public:
	using internals::PrimitiveInternalCoordinates::PrimitiveInternalCoordinates;
	MOCK_CONST_METHOD1(calc, scon::mathmatrix<coords::float_type>(coords::Representation_3D const& xyz));
	MOCK_CONST_METHOD2(calc_diff, scon::mathmatrix<coords::float_type>(coords::Representation_3D const& lhs, coords::Representation_3D const& rhs));
	MOCK_CONST_METHOD1(guess_hessian, scon::mathmatrix<coords::float_type>(internals::CartesianType const&));
	MOCK_METHOD1(Bmat, scon::mathmatrix<coords::float_type>& (internals::CartesianType const& cartesians));
	MOCK_METHOD1(Gmat, scon::mathmatrix<coords::float_type>& (internals::CartesianType const& cartesians));
	MOCK_METHOD1(transposeOfBmat, scon::mathmatrix<coords::float_type>(internals::CartesianType const& cartesians));
	MOCK_METHOD1(pseudoInverseOfGmat, scon::mathmatrix<coords::float_type>(internals::CartesianType const& cartesians));
};

class PrimitiveInternalSetTest : public testing::Test {
public:
	PrimitiveInternalSetTest();
	InternalCoordinates::CartesiansForInternalCoordinates cartesians;
	std::shared_ptr<internals::PrimitiveInternalCoordinates> testSystem;
	std::vector<std::vector<std::size_t>> molecules;
	ic_util::Graph<ic_util::Node> systemGraph;

	void createDistances();
	void createAngles();
	void createDihedrals();
	void createTranslations();
	void createRotations();

	void distanceCreationTest();
	void bondAngleCreationTest();
	void dihedralCreationTest();
	void tarnslationXCreationTest();
	/* void tarnslationYCreationTest();
	 void tarnslationZCreationTest();*/
	void rotationsCreationTest();
};

class MatricesTest : public testing::Test {
public:
	MatricesTest();
	InternalCoordinates::CartesiansForInternalCoordinates cartesians;
	std::shared_ptr<internals::PrimitiveInternalCoordinates> testSystem;
	std::vector<std::vector<std::size_t>> molecules;
	ic_util::Graph<ic_util::Node> systemGraph;

	void bMatrixTest();
	void gMatrixTest();
	void hessianGuessTest();
	void calculatePrimitiveInternalValuesTest();
};

class DelocalizedMatricesTest : public testing::Test {
public:
	DelocalizedMatricesTest();
	void delocalizedMatrixTest();
	void delocalizedBMatrixTest();
	void delocalizedGMatrixTest();
	void delocalizedInitialHessianTest();
	void internalDifferencesTest();
	void internalValuesForTricTest();


	InternalCoordinates::CartesiansForInternalCoordinates cartesians;
	std::shared_ptr<internals::TRIC> testSystem;
	std::vector<std::vector<std::size_t>> molecules;
	ic_util::Graph<ic_util::Node> systemGraph;
};

class ConverterMatricesTest : public testing::Test {
public:
	ConverterMatricesTest();
	void calculateInternalGradsTest();
	void applyInternalChangeTest();
	InternalCoordinates::CartesiansForInternalCoordinates cartesians;
	MockPrimitiveInternals testSystem;
	internals::InternalToCartesianConverter converter;
};

#endif

#endif