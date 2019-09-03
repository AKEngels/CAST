#include "BondGraph.h"

namespace ic_util{

	BondGraph::BondGraph(std::vector<std::string> const& symbols, coords::Representation_3D const& coordinates) : graph{ create_graph(symbols, coordinates) } {}

	std::vector<AtomNode> BondGraph::makeAtomNodeVector(std::vector<std::string> const& symbols) {
		std::vector<AtomNode> result;
		std::size_t enumerator = 0u;
		std::transform(symbols.cbegin(), symbols.cend(), std::back_inserter(result), [&enumerator](auto r) {
			return AtomNode{ ++enumerator, std::move(r) };
		});
		return result;
	}

	BondGraph::GraphType BondGraph::getGraph(std::vector<AtomNode> nodes, std::vector<std::pair<std::size_t, std::size_t>> connectivity) {
		GraphType result(connectivity.cbegin(), connectivity.cend(), nodes.size());

		for (auto i = 0u; i < nodes.size(); ++i) {
			result[i] = std::move(nodes[i]);
		}

		return result;
	}

	BondGraph::GraphType BondGraph::create_graph(std::vector<std::string> const& symbols, coords::Representation_3D const& coordinates) {
		return getGraph(makeAtomNodeVector(symbols), getConnectivity(symbols, coordinates));
	}

	std::vector<std::pair<std::size_t, std::size_t>> BondGraph::getConnectivity(std::vector<std::string> const& elem_vec, coords::Representation_3D const& cp_vec) {
		AtomConnector atomCreator(elem_vec, cp_vec);
		return atomCreator();
	}

	BondGraph::DistanceCreator::DistanceCreator(GraphType const& parentGraph) :
		source{ 0u },
		target{ 0u },
		edgeIterators{ boost::edges(parentGraph) },
		graph{ parentGraph }
	{}

	bool BondGraph::DistanceCreator::nextEdgeDistances() {
		if (edgeIterators.first == edgeIterators.second) return false;
		source = boost::source(*edgeIterators.first, graph);
		target = boost::target(*edgeIterators.first, graph);
		++edgeIterators.first;
		return true;
	}

	std::vector<std::tuple<std::size_t, std::size_t>> BondGraph::DistanceCreator::getBonds() {
		std::vector<std::tuple<std::size_t, std::size_t>> result;
		while (nextEdgeDistances()) {
			result.emplace_back( std::make_tuple( graph[source], graph[target] ));
		}
		return result;
	}

	BondGraph::AngleCreator::AngleCreator(GraphType const& parentGraph) :
		leftAtom{ 0u },
		middleAtom{ 0u },
		rightAtom{ 0u },
		vertexIterators{ boost::vertices(parentGraph) },
		resultPtr(nullptr),
		graph(parentGraph)
	{}

	std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> BondGraph::AngleCreator::getAngles() {
		std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> result;
		resultPtr = &result;
		while (nextVertex()) {
			addAngleForAllNeighbors();
		}
		resultPtr = nullptr;
		return result;
	}

	bool BondGraph::AngleCreator::nextVertex() {
		if (vertexIterators.first == vertexIterators.second) return false;
		middleAtom = *vertexIterators.first;
		vertexIterators.first++;
		return true;
	}

	void BondGraph::AngleCreator::addAngleForAllNeighbors() {
		auto allNeighbors = boost::adjacent_vertices(middleAtom, graph);
		spanLeftAndRightNeighborsForAngle(allNeighbors);
	}

	void BondGraph::AngleCreator::spanLeftAndRightNeighborsForAngle(std::pair<GraphType::adjacency_iterator, GraphType::adjacency_iterator>& neighbors) {
		while (findLeftAtom(neighbors)) {
			auto copyOfLeftNeighbors = neighbors;
			while (findRightAtom(copyOfLeftNeighbors)) {
				resultPtr->emplace_back(std::make_tuple(graph[leftAtom], graph[middleAtom], graph[rightAtom]));

			}
		}
	}

	bool BondGraph::AngleCreator::findLeftAtom(std::pair<GraphType::adjacency_iterator, GraphType::adjacency_iterator>& neighbors) {
		if (neighbors.first == neighbors.second) return false;
		leftAtom = *neighbors.first;
		++neighbors.first;
		return true;
	}

	bool BondGraph::AngleCreator::findRightAtom(std::pair<GraphType::adjacency_iterator, GraphType::adjacency_iterator>& neighborsLeft) {
		if (neighborsLeft.first == neighborsLeft.second) return false;
		rightAtom = *neighborsLeft.first;
		++neighborsLeft.first;
		return true;
	}

	BondGraph::DihedralCreator::DihedralCreator(GraphType const& parentGraph) :
		DistanceCreator{ parentGraph },
		outerLeft{ 0u },
		outerRight{ 0u },
		resultPtr{ nullptr }
	{}

	std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> BondGraph::DihedralCreator::getDihedrals() {
		std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> result;
		resultPtr = &result;
		while (nextEdgeDistances()) {
			findLeftAndRightAtoms();
		}
		resultPtr = nullptr;
		return result;
	}

	void BondGraph::DihedralCreator::findLeftAndRightAtoms() {
		auto leftVertices = boost::adjacent_vertices(source, graph);
		while (findLeftAtoms(leftVertices)) {
			auto rightVertices = boost::adjacent_vertices(target, graph);
			while (findRightAtoms(rightVertices)) {
				resultPtr->emplace_back(std::make_tuple(graph[outerLeft], graph[source], graph[target], graph[outerRight]));
			}
		}
	}

	bool BondGraph::DihedralCreator::findLeftAtoms(std::pair<GraphType::adjacency_iterator, GraphType::adjacency_iterator>& sourceNeighbors) {
		if (sourceNeighbors.first == sourceNeighbors.second) return false;
		outerLeft = *sourceNeighbors.first;
		++sourceNeighbors.first;
		if (outerLeft == target) return findLeftAtoms(sourceNeighbors);
		return true;
	}

	bool BondGraph::DihedralCreator::findRightAtoms(std::pair<GraphType::adjacency_iterator, GraphType::adjacency_iterator>& targetNeighbors) {
		if (targetNeighbors.first == targetNeighbors.second) return false;
		outerRight = *targetNeighbors.first;
		++targetNeighbors.first;
		if (outerRight == source) return findRightAtoms(targetNeighbors);
		return true;
	}

}