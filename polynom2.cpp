#include "polynom2.h"
#include <iostream>

polynom2::polynom2(int seed, bool debug) : mt_(seed), distribution_(-1.0, 1.0), debug_(debug)
{
}

complex_t 
polynom2::get_random_point()
{
	double real = distribution_(mt_);
	double imag = distribution_(mt_);
	return complex_t(real, imag);
}

complex_t 
polynom2::get_random_unit()
{
	auto theta = distribution_(mt_) * 3.14159;
	return complex_t(cos(theta), sin(theta));
}

double
polynom2::find_roots(int size)
{
	polynom_t p0(size);
	for (int i = 0; i < p0.size() - 1; i++)
		p0[i] = get_random_point();
	p0[p0.size() - 1] = 1; // leading coeff = 1

	std::vector<complex_t> roots;

	polynom_t q = p0;
	long long n = p0.size();
	do
	{
		complex_t z;
		// Newton solve
		complex_t delta;
		complex_t v;
		complex_t d;

		int max_try_start = 100;
		do
		{
			int max_steps = 1000;
			// Starting point
			z = get_random_unit();
			do
			{
				if (debug_)
					std::cout << z << std::endl;
				// Eval v = p(z) & d = p'(z)
				v = q[n - 1];
				d = 0;
				for (long long i = n - 2; i >= 0; i--)
				{
					v = v * z + q[i];
					d = d * z + ((double)(i + 1)) * q[i + 1];
				}
				delta = v / d;
				z = z - delta;
			} while (is_above(delta, 1e-16) && --max_steps>0); // check condition

			if (debug_ && max_steps <= 0)
			{
				std::cout << "==max step reached" << std::endl;
			}

		} while (is_above(v, 1e-6) && --max_try_start>0);

		if (debug_ && max_try_start <= 0)
			std::cout << "==== Solution approx" << std::endl;

		/*auto l2 = N2(eval(p0, z));
		std::cout << l2 << std::endl;
		if (l2 >= 1)
		{
			int a = 1;
		}*/

		roots.push_back(z);




		//  q = q/(x-z)
		v = q[n - 1];
		q[n - 1] = 0;
		for (long long i = n - 2; i >= 0; i--)
		{
			auto c = q[i];
			q[i] = v;
			v = v * z + c;
		}
		n--;
		if (debug_)
			multiply(q, n, z);
	} while (n >= 2);

	return distance_product(roots, p0);
}

complex_t 
polynom2::eval(const polynom_t& p, complex_t z)
{
	complex_t v = p[p.size() - 1];
	for (long long i = p.size() - 2; i >= 0; i--)
		v = v * z + p[i];
	return v;
}

double 
polynom2::N2(const complex_t& z)
{
	return z.real() * z.real() + z.imag() * z.imag();
}

bool 
polynom2::is_above(const complex_t& z, double epsilon)
{
	return N2(z) > epsilon;
}


void 
polynom2::multiply(const polynom_t& q, long long n, complex_t c)
{
	std::cout << 0 << ": " << -q[0] * c << std::endl;
	for (long long i = 1; i < n; i++)
		std::cout << i << ": " << q[i - 1] - q[i] * c << std::endl;
	std::cout << n << ": " << q[n - 1] << std::endl;
}

// Regenerate polynom: product(X-ri)
double
polynom2::distance_product(const std::vector<complex_t>& roots, const polynom_t& p0)
{
	polynom_t q(roots.size() + 1);
	q[0] = -roots[0];
	q[1] = 1;
	// q(z) <= q(z).(z-ri)
	for (size_t j = 1; j < roots.size(); j++)
	{
		q[j + 1] = q[j];
		for (long long i = j; i >= 1; i--)
		{
			q[i] = q[i - 1] - q[i] * roots[j];
		}
		q[0] = -q[0] * roots[j];
	}
	if (p0.size() != q.size())
		throw std::runtime_error("Sizes are differents: " + std::to_string(p0.size()) + " " + std::to_string(q.size()));
	double diff = 0;
	for (size_t j = 0; j < p0.size(); j++)
		diff += N2(p0[j] - q[j]);
	return diff;
}
