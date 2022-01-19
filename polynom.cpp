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

// Show p(z).(z-c)
void multiply(const polynom_t& q, long long n, complex_t c)
{
	std::cout << 0 << ": " << -q[0] * c << std::endl;
	for (long long i = 1; i < n ; i++)
		std::cout << i << ": " << q[i - 1] - q[i] * c << std::endl;
	std::cout << n << ": " << q[n - 1] << std::endl;
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
	complex_t z0(0.0, 0.0);

	polynom_t q = p0;
	long long n = p0.size();
	do
	{
		auto z = z0;
		// Newton solve
		complex_t delta;
		do
		{
			std::cout << z << std::endl;
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
		} while (is_above(delta, 1e-16)); // check condition
		std::cout << "root: " << z << std::endl;

		//  q = q/(x-z)
		complex_t v = q[n - 1];
		q[n - 1] = 0;
		for (long long i = n - 2; i >= 0; i--)
		{
			auto c = q[i];
			q[i] = v;
			v = v * z + c;
		}
		n--;
		multiply(q, n, z);
	} 
	while (n >= 2);

	return 1234.0;

}