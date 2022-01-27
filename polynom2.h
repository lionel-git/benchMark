#pragma once

#include <random>
#include <complex>

typedef std::complex<double> complex_t;
typedef std::vector<std::complex<double>> polynom_t;

class polynom2
{
public:
	polynom2(int seed, bool debug = false);
	double find_roots(int n);

	// Static utils
private:
	static complex_t eval(const polynom_t& p, complex_t z);
	static double N2(const complex_t& z);
	static bool is_above(const complex_t& z, double epsilon);
	static void multiply(const polynom_t& q, long long n, complex_t c);
	static double distance_product(const std::vector<complex_t>& roots, const polynom_t& p0);

	complex_t get_random_point();

private:
	std::mt19937 mt_;
	std::uniform_real_distribution<double> distribution_;
	bool debug_;
};