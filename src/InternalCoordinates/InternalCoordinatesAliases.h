#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATESALISES_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATESALISES_H_

#include<memory>

namespace ic_util {
	enum class period;
	class BondGraph;
}

namespace internals {
	using float_type = double;
}

namespace scon {/*
	template<typename T>
	class c3;*/
	template <typename T> 
	class mathmatrix;
	/*template<typename T, typename Allocator>
	class vector;*/
}
/*
namespace coords {
	using r3 = scon::c3<internals::float_type>;
	using Cartesian_Point = r3;
	template<typename T, typename Allocator> using Container = scon::vector<T, Allocator>;
	
	//using Representation_3D = Container<Cartesian_Point, std::allocator<internals::float_type>>;

	class Coordinates;
}
*/

namespace InternalCoordinates {
	struct InternalCoordinate;
	class Rotator;
	class temporaryCartesian;
	class CartesiansForInternalCoordinates;
	class InternalCoordinatesBuilder;
}

namespace internals {
	using CartesianType = InternalCoordinates::CartesiansForInternalCoordinates;
	class AppropriateStepFinder;
	class InternalToCartesianConverter;
}

#endif // CAST_INTERNALCOORDINATES_INTERNALCOORDINATESALISES_H_