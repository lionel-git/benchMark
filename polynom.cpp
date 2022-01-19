#include <vector>
#include <complex>
#include <random>
#include <iostream>

typedef std::complex<double> complex_t;
typedef std::vector<std::complex<double>> polynom_t;

// return P(z)
static complex_t eval(const polynom_t& p, complex_t z)
{
	complex_t v = p[p.size() - 1];
	for (long long i = p.size() - 2; i >= 0; i--)
		v = v * z + p[i];
	return v;
}

bool is_above(complex_t z, double epsilon)
{
	return z.real() * z.real() + z.imag() * z.imag() > epsilon;
}

double find_roots(int size)
{
	std::mt19937 mt(1234);
	std::uniform_real_distribution<double> distribution(-1.0, 1.0);

	polynom_t p0(size);
	for (int i = 0; i < p0.size() - 1; i++)
	{
		p0[i].real(distribution(mt));
		p0[i].imag(distribution(mt));
	}
	p0[p0.size() - 1] = 1; // leading coeff = 1

	// Starting point
	complex_t z(distribution(mt), distribution(mt));

	polynom_t q = p0;
	long long n = p0.size();

	// Newton solve
	complex_t delta;
	do
	{
		// Eval v = p(z) & d = p'(z)
		complex_t v = q[n - 1];
		complex_t d = 0;
		for (long long i = n - 2; i >= 0; i--)
		{
			v = v * z + q[i];
			d = d * z + ((double)(i + 1)) * q[i + 1];
		}
		delta = v / d;
		z = z - delta;
		std::cout << z << std::endl;
	}
	while (is_above(delta, 1e-16)); // check condition
	
	//  p = p/(x-z)
	complex_t v = q[n - 1];
	q[n - 1] = 0;
	for (long long i = n - 2; i >= 0; i--)
	{
		q[i] = v;
		v = v * z + q[i];
	}	
	return 1234.0;

}