#ifndef CAST_INTERNALCOORDINTATES_BONDGRAPH_BONDGRAPH_H_
#define CAST_INTERNALCOORDINTATES_BONDGRAPH_BONDGRAPH_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/variant/get.hpp>
#include <boost/graph/graphviz.hpp>
#include <string>
#include <utility>
#include <vector>

#include "../../coords.h"
#include "../../coords_rep.h"

#include "AtomConnector.h"

/*!
 *  \addtogroup ic_util
 *  @{
 */
namespace ic_util {

	/*!
	\details The Node struct and Graph class can be thought of as instantiations of
	the graph abstraction provided by the Boost Graph Library (BGL).
	\see "The Boost Graph Library", Jeremy Siek et al., Addison-Wesley, 2002
	*/

	/*!
	\brief Vertex class.
	\details User-defined class type for the VertexProperty template argument of
	boost::adjacency_list<>. The members of Node are the properties (a.k.a. labels)
	that each vertex possesses; in BGL terminology this is called bundled property.
	*/
	struct AtomNode {
		unsigned int atom_serial;
		std::string element;
	};

	/*!
	\brief User-defined Graph class that represents an actual instantiation of BGL's
	graph abstraction.
	\tparam Atom_type User-defined type that encapsulates atom-specific data
	retrieved by the PDB parser.
	*/
	class BondGraph {

	public:
		BondGraph(std::vector<std::string> const& symbols, coords::Representation_3D const& coordinates);

		/*!
		\brief User-defined graph type.
		\details Basis of the graph is either an adjacency list or a matrix.
		\tparam vecS The first template param is the container used to store the edge
		list; in this case std::vector.
		\tparam vecS The second template param is the container used to store the
		vertex list; again std::vector is used.
		\tparam undirectedS Specifies that the graph is undirected.
		\tparam Node Specifies vertex properties.
		*/
		using GraphType = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, AtomNode>;
	private:
		GraphType graph; 

		
		static std::vector<AtomNode> makeAtomNodeVector(std::vector<std::string> const& symbols);
		static GraphType getGraph(std::vector<AtomNode> nodes, std::vector<std::pair<std::size_t, std::size_t>> connectivity);
		static GraphType create_graph(std::vector<std::string> const& symbols, coords::Representation_3D const& coordinates);
		static std::vector<std::pair<std::size_t, std::size_t>> getConnectivity(std::vector<std::string> const& elem_vec, coords::Representation_3D const& cp_vec);

		class DistanceCreator {
		public:
			DistanceCreator(GraphType const& graph);
			virtual ~DistanceCreator() = default;

			std::vector<std::tuple<std::size_t, std::size_t>> getBonds();

		protected:
			bool nextEdgeDistances();
			std::size_t source, target;
			std::pair<GraphType::edge_iterator, GraphType::edge_iterator> edgeIterators;
			GraphType const& graph;
		};

		class AngleCreator {
		public:
			AngleCreator(GraphType const& graph);
			virtual ~AngleCreator() = default;

			std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> getAngles();

		protected:
			bool nextVertex();
			void addAngleForAllNeighbors();
			void spanLeftAndRightNeighborsForAngle(std::pair<GraphType::adjacency_iterator, GraphType::adjacency_iterator>& neighbors);
			bool findLeftAtom(std::pair<GraphType::adjacency_iterator, GraphType::adjacency_iterator>& neighbors);

			bool findRightAtom(std::pair<GraphType::adjacency_iterator, GraphType::adjacency_iterator>& neighborsLeft);

			std::size_t leftAtom, middleAtom, rightAtom;
			std::pair<GraphType::vertex_iterator, GraphType::vertex_iterator> vertexIterators;
			std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> * resultPtr;
			GraphType const& graph;
		};

		class DihedralCreator : public DistanceCreator {
		public:
			DihedralCreator(GraphType const& graph);
			virtual ~DihedralCreator() = default;

			std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> getDihedrals();
		protected:
			void findLeftAndRightAtoms();
			bool findLeftAtoms(std::pair<GraphType::adjacency_iterator, GraphType::adjacency_iterator>& sourceNeighbors);
			bool findRightAtoms(std::pair<GraphType::adjacency_iterator, GraphType::adjacency_iterator>& targetNeighbors);

			std::size_t outerLeft, outerRight; 
			std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> * resultPtr;

		};

		std::vector<std::vector<std::size_t>> molecules;

		void findMolecules() {
			std::vector<std::size_t> components(boost::num_vertices(graph));
			auto numberOfMolecules = boost::connected_components(graph, components.data());
			molecules.resize(numberOfMolecules);
			for (auto i = 0u; i < components.size(); ++i) {
				molecules[components[i]].emplace_back(graph[i].atom_serial);
			}
		}
		

	public:

		AtomNode const& getAtom(std::size_t const i) const { return graph[i]; }
		std::size_t getNumberOfAtoms() const { return boost::num_vertices(graph); }
		
		std::vector<std::vector<std::size_t>> const& getMolecules() {
			if(molecules.empty()) findMolecules();
			return molecules;
		}
		
		void resetMolecules() {molecules.clear();}

		std::vector<std::tuple<std::size_t, std::size_t>> getBonds() const { return DistanceCreator(graph).getBonds(); }

		std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> getAngles() const { return AngleCreator(graph).getAngles(); }

		std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> getDihedrals() const { return DihedralCreator(graph).getDihedrals(); }


		/*!
		\brief Creates the graph data structure. Called during instantiation of the
		Graph class.
		\param b_at Vector of bonded atom pairs.
		\param vec Vector of atoms.
		\return Graph_type object.
		*/
		

		/*!
		\brief Prints the edges of a graph to std::cout.
		\param graph Graph, whose edges shall be printed.
		\return The two vertices of an edge are printed.
		*/
		std::ostream & printEdges(std::ostream & os) {
			using boost::edges;
			using boost::source;
			using boost::target;

			os << "Bond graph edges: \n";
			auto ed = edges(graph);
			for (auto it = ed.first; it != ed.second; ++it) {
				auto u = source(*it, graph);
				auto v = target(*it, graph);
				os << graph[u].atom_serial << " -> " << graph[v].atom_serial << "\n";
			}
		}

		/*!
		\brief Prints the vertices of a graph to std::cout.
		\param graph Graph, whose vertices shall be printed.
		\return Each vertex with its associated properties is printed.
		*/
		std::ostream & print_vertices(std::ostream & os) {
			using boost::vertices;

			os << "Bond graph vertices: \n";
			auto vert = vertices(graph);
			for (auto it = vert.first; it != vert.second; ++it) {
				os << "Atom: " << graph[*it].atom_serial << "; " << graph[*it].element << "\n";
			}
		}
	};
}

/*! @} End of ic_util group*/

#endif // CAST_INTERNALCOORDINTATES_BONDGRAPH_H_
