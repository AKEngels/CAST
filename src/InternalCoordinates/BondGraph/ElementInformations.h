#ifndef CAST_INTERNALCOORDINATES_BONDGRAPH_ELEMNTINFORMATION_H_
#define CAST_INTERNALCOORDINATES_BONDGRAPH_ELEMNTINFORMATION_H_

#include<string>

#include "../InternalCoordinatesAliases.h"

namespace ic_util{

enum class period : int {
	  none = 0,
	  one,
	  two,
	  three
  };

std::size_t element_number(std::string const& key);

internals::float_type element_radius(std::string const& key);

period element_period(std::string const& key);

}

#endif // CAST_INTERNALCOORDINATES_BONDGRAPH_ELEMNTINFORMATION_H_