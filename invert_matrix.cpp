#include <memory>
#include <iostream>
#include <random>
#include <iomanip>

static bool debug_matrix = false;

typedef std::vector<std::vector<double>> matrix_t;

static double spread_identity_multiply_matrix(const matrix_t& m1, int shift_colm1, const matrix_t& m2, int shift_colm2)
{
	double distance = 0.0;
	for (size_t i = 0; i < m1.size(); i++)
	{
		for (size_t j = 0; j < m1.size(); j++)
		{
			double s = 0.0;
			for (size_t k = 0; k < m1.size(); k++)
				s += m1[i][k + shift_colm1] * m2[k][j + shift_colm2];
			if (i == j)
				distance += fabs(s - 1.0);
			else
				distance += fabs(s);
		}
	}
	return distance;
}

static void show_matrix(const matrix_t& m)
{
	if (debug_matrix)
	{
		std::cout << "=========================" << std::endl;
		for (size_t i = 0; i < m.size(); i++)
		{
			for (size_t j = 0; j < 2 * m.size(); j++)
			{
				std::cout << std::left << std::setprecision(3) << std::setw(8) << m[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

double invert_matrix(size_t size)
{
	std::mt19937 mt(1234);
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	matrix_t m;
	for (size_t i = 0; i < size; i++)
	{
		std::vector<double> row(2 * size);
		for (size_t j = 0; j < size; j++)
		{
			row[j] = distribution(mt);
			if (i == j)
				row[j + size] = 1.0;
			else
				row[j + size] = 0.0;
		}
		m.push_back(row);
	}
	show_matrix(m);

	auto m0 = m;

	size_t pivot = 0;

	while (pivot < size)
	{
		// Search max abs value on pivot column
		double max = fabs(m[pivot][pivot]);
		size_t row_max = pivot;
		for (size_t i = pivot + 1; i < size; i++)
		{
			if (fabs(m[i][pivot]) > max)
			{
				max = fabs(m[i][pivot]);
				row_max = i;
			}
		}
		if (max < 1e-8)
			throw std::runtime_error("Matrix is not invertible!");

		if (row_max != pivot)
			std::swap(m[row_max], m[pivot]);
		show_matrix(m);

		double v = m[pivot][pivot];
		for (size_t i = 0; i < size; i++)
		{
			if (i != pivot)
			{
				double mult = m[i][pivot] / v;
				for (size_t j = 0; j < 2 * size; j++)
					m[i][j] -= m[pivot][j] * mult;
			}
		}
		show_matrix(m);
		pivot++;
	}

	for (size_t i = 0; i < size; i++)
	{
		double mult = 1.0 / m[i][i];
		for (size_t j = 0; j < 2 * size; j++)
			m[i][j] *= mult;
	}
	show_matrix(m);

	show_matrix(m0);
	show_matrix(m);
	return spread_identity_multiply_matrix(m0, 0, m, size);
}
