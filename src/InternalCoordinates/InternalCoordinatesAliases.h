#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATESALISES_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATESALISES_H_

namespace ic_util {
	enum class period;
	class BondGraph;
}

namespace internals {
	using float_type = double;
}

namespace scon {
	template<typename T>
	class c3;
	template <typename T> 
	class mathmatrix;
}

namespace coords {
	using r3 = scon::c3<internals::float_type>;
	using Cartesian_Point = r3;
}


namespace InternalCoordinates {
	struct InternalCoordinate;
	class Rotator;
	class temporaryCartesian;
	class CartesiansForInternalCoordinates;
	class InternalCoordinatesBuilder;
}

#endif // CAST_INTERNALCOORDINATES_INTERNALCOORDINATESALISES_H_