/**
CAST 3
md.h
Purpose: header for molecular dynamics simulation

@version 1.0
*/

#pragma once 

#include <vector>
#include <string>



/**
*namespace for debug and dev stuff that is part of the MD Implementation but of no concern for routine usage
*/
namespace md
{

	/**struct that contains information about an atom pair that is to be analyzed*/
	struct ana_pair
	{
		/**index of first atom (starting with 0)*/
		int a;
		/**index of second atom (starting with 0)*/
		int b;
		/**element symbol of first atom*/
		std::string symbol_a;
		/**element symbol of second atom*/
		std::string symbol_b;
		/**name of first atom (element and tinker atom index)*/
		std::string name_a;
		/**name of second atom (element and tinker atom index)*/
		std::string name_b;
		/**legend of this atom pair in the graph*/
		std::string legend;
		/**distances for every MD frame*/
		std::vector<double> dists;

		/**constructor
		@param p1: tinker atom index of first atom (i.e. starting with 1)
		@param p2: tinker atom index of second atom (i.e. starting with 1)*/
		ana_pair(int p1, int p2) { a = p1 - 1; b = p2 - 1; }

		/**returns a string with information about the atom pair*/
		std::string info()
		{
			std::string result = "Atoms: " + name_a + " , " + name_b + "\n";
			for (auto d : dists)
			{
				result += std::to_string(d) + " , ";
			}
			return result + "\n";
		}
	};

	/**information about a zone for which temperature is to by analyzed*/
	struct zone
	{
		/**legend for plotting*/
		std::string legend;
		/**atom indizes (starting with 0)*/
		std::vector<int> atoms;
		/**temperatures for every MD step*/
		std::vector<double> temperatures;
	};

}

