#include "interpolation.h"
#include <iostream>


Lagrange_interp::Lagrange_interp(std::vector<double> const& gridx, std::vector<double> const& gridy) : Base_interpolation(gridx, gridy),
pointsx(gridx), pointsy(gridy)
{
	grid_size = gridx.size();
	lambda.resize(grid_size);
	for (size_t i = 0; i < lambda.size(); i++) lambda[i].resize(grid_size);
	build();
}

void Lagrange_interp::build()
{
	lambda[0][0] = 1.0;
	size_t k(1);
	while (k < grid_size)
	{
		for (size_t i = 0; i < k; i++)
		{
			lambda[k][i] = lambda_call(k, i);
			lambda[k][k] -= lambda_call(k, i);
		}
		k++;
	}
}

double Lagrange_interp::lambda_call(size_t k, size_t i)
{
	return lambda[k - 1][i] / (pointsx[i] - pointsx[k]);
}

double Lagrange_interp::interpolate(double step)
{
	double sum(0.0), sum2(0.0);
	for (size_t i = 0; i < lambda[grid_size - 1].size(); i++)
	{
		if (step == pointsx[i]) return pointsy[i];
		sum += (lambda[grid_size - 1][i] / (step - pointsx[i])) * pointsy[i];
		sum2 += lambda[grid_size - 1][i] / (step - pointsx[i]);
	}
	return (sum / sum2);
}

Linear_interp_sorted::Linear_interp_sorted(std::vector<double> const& gridx, std::vector<double> const& gridy) : Base_interpolation(gridx, gridy),
pointsx(gridx), pointsy(gridy)
{
	grid_size = gridx.size();
}

double Linear_interp_sorted::interpolate(double step)
{
	size_t holder(Linear_interp_sorted::Base_interpolation::locate(step, pointsx));
	return (pointsy[holder] + pointsy[holder + 1]) / 2;
}



Spline_interp_natural::Spline_interp_natural(std::vector<double> const& gridx, std::vector<double> const& gridy) : Base_interpolation(gridx, gridy),
pointsx(gridx), pointsy(gridy)
{
	grid_size = gridx.size();
	div2.resize(gridx.size());
}

size_t Base_interpolation::locate(double x, std::vector<double>& gridx)
{
	size_t pos(0), upper(grid_size - 1), lower(0), midpoint(0);
	if (gridx[0] <= gridx[grid_size - 1]) asc = true;
	else asc = false;
	while ((upper - lower) > 1)
	{
		midpoint = (upper + lower) / 2;
		if ((x >= gridx[midpoint]) == asc) lower = midpoint;
		else upper = midpoint;
	}
	if (x == gridx[0])pos = 0;
	else if (x == gridx[grid_size]) pos = grid_size - 1;
	else pos = lower;
	return pos;

}

double Spline_interp_natural::interpolate(double step)
{
	double sigma(0.0), decomp_p(0.0), upper_bound_q(0.0), upper_bound_n(0.0), h(0.0),
		a(0.0), b(0.0), y(0.0);
	std::vector<double>decompo_u(grid_size, 0.0);
	div2[0] = 0.0;
	for (size_t i = 1; i < grid_size - 1; i++)
	{
		sigma = (pointsx[i] - pointsx[i - 1]) / (pointsx[i + 1] - pointsx[i - 1]);
		decomp_p = sigma * div2[i - 1] + 2.0;
		div2[i] = (sigma - 1.0) / decomp_p;
		decompo_u[i] = (6.0 * ((pointsy[i + 1] - pointsy[i]) / (pointsx[i + 1] - pointsx[i]) - (pointsy[i] - pointsy[i - 1])
			/ (pointsx[i] - pointsx[i - 1])) / (pointsx[i + 1] - pointsx[i - 1]) - sigma * decompo_u[i - 1]) / decomp_p;

	}
	div2[grid_size - 1] = (upper_bound_n - upper_bound_q * decompo_u[grid_size - 2]) / (upper_bound_n * div2[grid_size - 2] + 1.0);
	for (size_t l = grid_size - 2; l > 0; l--) div2[l] = div2[l] * div2[l + 1] + decompo_u[l];

	size_t lower(Spline_interp_natural::Base_interpolation::locate(step, pointsx));
	h = pointsx[lower + 1] - pointsx[lower];
	if (h == 0) throw std::runtime_error("bad input to cubic-spline routine");
	a = (pointsx[lower + 1] - step) / h;
	b = (step - pointsx[lower]) / h;
	y = a * pointsx[lower] + b * pointsy[lower + 1] + ((a * a * a - a) * div2[lower] + (b * b * b - b) * div2[lower + 1]) * (h * h) / 6.0;
	return y;
}
