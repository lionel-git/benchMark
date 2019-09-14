#include <iostream>
#include <chrono>

long iterate(double cx, double cy, int max)
{
	int k = 0;
	double x = 0.0;
	double y = 0.0;
	double x2 = x * x;
	double y2 = y * y;
	while (k < max && (x2 + y2 < 4.0))
	{
		// z = z^2 + c
		y = 2 * x * y + cy;
		x = x2 - y2 + cx;
		x2 = x * x;
		y2 = y * y;
		k++;
	}
	return k;
}


long benchMandel(double dx, double dy)
{
	long total = 0;
	int max = 200;
	for (double x = -2.0; x < 2.0; x += dx)
		for (double y = -2.0; y < 2.0; y += dy)
		{
			total += iterate(x, y, max);
		}
	return total;
}

// Avoid optimiser to call start/end before benchMark runs ...
#ifdef __linux__
std::chrono::time_point<std::chrono::_V2::system_clock> get_time(long b)
#else
std::chrono::time_point<std::chrono::steady_clock> get_time(long b)
#endif
{
	if (b == 1)
		return std::chrono::high_resolution_clock::now();
	else
		return std::chrono::high_resolution_clock::now();
}

int main(int argc, char** argv)
{  
	auto start = get_time(0);
	long total = benchMandel(0.0005, 0.0005);
	auto end = get_time(total);
	std::cout << "res=" << total << std::endl;
	std::chrono::duration<double> diff = end - start;
	std::cout << diff.count() << " s" << std::endl;
}
