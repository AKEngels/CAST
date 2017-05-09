#ifndef H_2DSCAN
#define H_2DSCAN

#include<vector>
#include<sstream>
#include<iterator>
#include<memory>

#include"coords.h"
#include"coords_rep.h"
#include"configuration.h"

class Scan2D {
public:
	Scan2D(coords::Coordinates & coords);

	


private:

	struct what {
		std::string what_kind;
		std::vector<std::size_t> atoms;
		double change;
		double to_position;
		std::size_t scans;
	};
	class Input_types {
	public:
		virtual void fill_what(std::vector<std::string> & splitted_vals, Scan2D::what & to_fill) = 0;
	};
	class Normal_Input : public Input_types {
		void fill_what(std::vector<std::string> & splitted_vals, Scan2D::what & to_fill) override;
	};

	void parse_input(std::vector<std::string> const & input);
	Input_types * Input_Factory(std::size_t size, std::string kind);

	what x_axis;
	what y_axis;

	coords::Coordinates & _coords;
};

#endif