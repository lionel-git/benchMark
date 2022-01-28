#pragma once

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <processthreadsapi.h>

bool check_avx()
{
	return IsProcessorFeaturePresent(PF_AVX_INSTRUCTIONS_AVAILABLE);
}

bool check_avx2()
{
	return IsProcessorFeaturePresent(PF_AVX2_INSTRUCTIONS_AVAILABLE);
}

bool check_avx512()
{
	return IsProcessorFeaturePresent(PF_AVX512F_INSTRUCTIONS_AVAILABLE);
}

#elif defined(__arm__)

bool check_avx()
{
	return false;
}

bool check_avx2()
{
  	return false;
}

bool check_avx512()
{
  	return false;
}

#else

// https://gcc.gnu.org/onlinedocs/gcc/x86-Built-in-Functions.html

bool check_avx()
{
	return __builtin_cpu_supports("avx");
}

bool check_avx2()
{
	return __builtin_cpu_supports("avx2");
}

bool check_avx512()
{
	return __builtin_cpu_supports("avx512f");
}

#endif


