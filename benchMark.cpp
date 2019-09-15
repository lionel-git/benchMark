#include <iostream>
#include <chrono>
#include <cstring>
#include <math.h>

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

long benchHeat(double dx, double dy)
{
	int nx = (int)10.0 / dx;
	int ny = (int)10.0 / dx;
	int L = nx * ny;
	if (L == 0)
		throw new std::runtime_error("L=0");
	std::cout << "Mem used: " << 2.0 * L * sizeof(double)/(1024.0*1024.0) << " Mo" << std::endl;
	double* values = new double[L];
	double* values2 = new double[L];

	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			values[i * ny + j] = 0.0;

	for (int j = 0; j < ny; j++)	
	{
		values[0 * ny + j] = j * 45.0;
		values[(nx-1) * ny + j] = - j * 4.0;
	}

	for (int i = 0; i < nx; i++)
	{
		values[i * ny + 0] = i * 17.0;
		values[i * ny + ny - 1] = -i * 5.0;
	}

	memcpy(values2, values, L*sizeof(double));

	double* values_src = values;
	double* values_dest = values2;

	double diff = 0.0;
	long k = 0;
	do
	{
		// Update temp
		diff = 0.0;
		for (int i = 1; i < nx - 1; i++)
			for (int j = 1; j < ny - 1; j++)
			{
				values_dest[i * ny + j] = 0.25 *
					(values_src[i * ny + (j - 1)] + values_src[i * ny + (j + 1)]
						+ values_src[(i - 1) * ny + j] + values_src[(i + 1) * ny + j]);
				diff += fabs(values_dest[i * ny + j] - values_src[i * ny + j]);
			}
		// switch vectors
		auto temp = values_dest;
		values_dest = values_src;
		values_src = temp;
		k++;
	} while (diff > 0.1*L);
	return k;
}

void bench(const std::string& func_name, double dx, double dy, long(*bench_func)(double ddx, double ddy))
{
	std::cout << "Start " << func_name << std::endl;
	auto start = get_time(0);
	long total = bench_func(dx, dy);
	auto end = get_time(total);
	std::cout << "res=" << total << std::endl;
	std::chrono::duration<double> diff = end - start;
	std::cout << diff.count() << " s" << std::endl;
	std::cout << "==============" << std::endl;
}

int main(int argc, char** argv)
{  
//	bench("benchMandel", 0.0005, 0.0005, benchMandel);
	bench("benchHeat", 0.01, 0.01, benchHeat);
} 
