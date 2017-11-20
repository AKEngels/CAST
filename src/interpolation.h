#include <vector>
#include <cstddef>


class Base_interpolation

{

protected:

	bool asc;
	size_t grid_size, order_k;
	
		Base_interpolation( std::vector<double> const & gridx, std::vector<double> const & gridy)
		{
			grid_size = gridx.size();
			if (gridx.size() != gridy.size()) throw ("Unequally sized grid");
			asc = false;
		};
	size_t locate(double x, std::vector<double> & gridx);
		

public:

	virtual double interpolate(double step) = 0;
	

};


class Lagrange_interp : Base_interpolation

{
public:

	Lagrange_interp(std::vector<double> const & gridx, std::vector<double> const & gridy);
	double interpolate(double step);

private:
	std::vector <double> pointsx, pointsy, weights;
	std::vector <std::vector <double> > lambda;
	void build();
	double lambda_call(size_t k, size_t i);
	
};


class Linear_interp_sorted : Base_interpolation

{

public:

	Linear_interp_sorted(std::vector<double> const & gridx, std::vector<double> const & gridy);
	double interpolate(double step);

private:

	std::vector <double> pointsx, pointsy;
	
};

class Spline_interp_natural : Base_interpolation
{
public:

	Spline_interp_natural(std::vector<double> const & gridx, std::vector<double> const & gridy);
	double interpolate(double step);


private:
	std::vector <double> pointsx, pointsy, div2 ;
	
};


