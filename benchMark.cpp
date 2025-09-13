#include <iostream>
#include <chrono>
#include <cstring>
#include <math.h>
#include <thread>
#include <vector>
#include <string>
#include <complex>
#include <random>

#include "features_check.h"

#include "polynom2.h"

bool G_pin_threads = true;

extern void FFT(int dir, long m, std::complex<double> x[]);
extern double invert_matrix(size_t size);

long iterate(double cx, double cy, int max)
{
    int k = 0;
    double x = 0.0;
    double y = 0.0;
    double x2 = x * x;
    double y2 = y * y;
    while ((k & 1) == 1 || (k < max && (x2 + y2 < 4.0)))
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

// unroll 2 iterations
// z = (z^2 + c)^2 + c = z^4 + 2.c.z^2 + c^2 + c
// z = z^4 + B.z^2 + C avec B=2c C=c^2+c
long iterate2(double cx, double cy, int max)
{
    int k = 0;
    double Cx = cx * cx - cy * cy + cx;
    double Cy = 2 * cx * cy + cy;
    double Bx = 2 * cx;
    double By = 2 * cy;
    double x = 0.0;
    double y = 0.0;
    double x2 = x * x;
    double y2 = y * y;
    while (k < max && (x2 + y2 < 4.0))
    {
        // Z2=z^2
        double Z2x = x2 - y2;
        double Z2y = 2 * x * y;

        // z = z^4 + B.z^2 + C = Z2.(Z2 + B) + C
//		x = (Z2x * Z2x - Z2y * Z2y) + (Bx * Z2x - By * Z2y) + Cx;
//		y = (2 * Z2x * Z2y) + Bx * Z2y + By * Z2x + Cy;

        // z= Z2 . K + C avec K=Z2 + B
        double Kx = Z2x + Bx;
        double Ky = Z2y + By;

        x = Z2x * Kx - Z2y * Ky + Cx;
        y = Z2x * Ky + Z2y * Kx + Cy;

        x2 = x * x;
        y2 = y * y;
        k += 2;
    }
    return k;
}

void benchMandel(double dx, double dy, long long& total)
{
    total = 0;
    int max = 200;
    for (double x = -2.0; x < 2.0; x += dx)
        for (double y = -2.0; y < 2.0; y += dy)
        {
            total += iterate(x, y, max);
        }
}

void benchMandel2(double dx, double dy, long long& total)
{
    total = 0;
    int max = 200;
    for (double x = -2.0; x < 2.0; x += dx)
        for (double y = -2.0; y < 2.0; y += dy)
        {
            total += iterate2(x, y, max);
        }
}


void benchHeat(double dx, double /* dy */, long long& ret)
{
    int nx = (int)(10.0 / dx);
    int ny = (int)(10.0 / dx);
    int L = nx * ny;
    if (L == 0)
        throw new std::runtime_error("L=0");
    std::cout << "Mem used: " << 2.0 * L * sizeof(double) / (1024.0 * 1024.0) << " Mo" << std::endl;
    double* values = new double[L];
    double* values2 = new double[L];

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            values[i * ny + j] = 0.0;

    for (int j = 0; j < ny; j++)
    {
        values[0 * ny + j] = j * 45.0;
        values[(nx - 1) * ny + j] = -j * 4.0;
    }

    for (int i = 0; i < nx; i++)
    {
        values[i * ny + 0] = i * 17.0;
        values[i * ny + ny - 1] = -i * 5.0;
    }

    memcpy(values2, values, L * sizeof(double));

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
    } while (diff > 0.1 * L);
    ret = k;
    delete[] values;
    delete[] values2;
}

void benchPi(double dx, double /* dy */, long long& ret)
{
    double sum = 0.0;
    for (double x = 0.0; x < 1.0; x += dx)
        sum += sqrt(1.0 - x * x);
    ret = (long long)(dx * sum * 40000000000.0);
}

void benchFft(double m0, double c, long long& total)
{
    long m = (long)(m0 + 0.5);
    if (m <= 1 || m >= 32)
    {
        std::cerr << "m is invalid: " << m << std::endl;
        return;
    }

    long N = 1l << m;
    std::unique_ptr<std::complex<double>[]> z(new std::complex<double>[N]);

    std::mt19937 mt(1234);

    for (long i = 0; i < N; i++)
    {
        auto v = mt();
        if (v & 1)
        {
            z[i].real(c);
            z[i].imag(1.0 - c);
        }
        else
        {
            z[i].real(1.0 - c);
            z[i].imag(c);
        }
    }

    FFT(+1, m, z.get());
    FFT(-1, m, z.get());

    double ret = 0.0;
    for (long i = 0; i < N; i++)
        ret += z[i].real() + z[i].imag();
    total = (long long)(10.0 * ret + 0.5);
}

void benchMatrix(double size, double a, long long& total)
{
    auto distance = invert_matrix((int)size);
    total = (long long)(distance * 1e10);
}

void benchPolynom(double size, double a, long long& total)
{
    double distance = 0.0;
    for (int i = 0; i < (int)a; i++)
    {
        polynom2 p(i);
        distance += p.find_roots((int)size);
    }
    total = (long long)(distance * 1e15);
}

// == Utils for bench ================================================

// Avoid optimiser to call start/end before benchMark runs ...
auto get_time(long long b)
{
    if (b == 1)
        return std::chrono::high_resolution_clock::now();
    else
        return std::chrono::high_resolution_clock::now();
}

void set_thread_affinity(std::thread::native_handle_type thread_handle, unsigned int core_id)
{
#if defined(_WIN32)
    const auto core_mask = 1ULL << core_id;
    const auto rc = SetThreadAffinityMask(thread_handle, core_mask);
    if (rc == 0)
        std::cerr << "Error calling SetThreadAffinityMask with core_mask=" << core_mask << ": " << GetLastError() << std::endl;
#else
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);
    int rc = pthread_setaffinity_np(thread_handle, sizeof(cpu_set_t), &cpuset);
    if (rc != 0) {
        std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
    }
    /*
    CPU_ZERO(&cpuset);
    int rc2 = pthread_getaffinity_np(thread_handle, sizeof(cpu_set_t), &cpuset);
    if (rc2 != 0) {
        std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
    }
    for (size_t i=0;i<CPU_SETSIZE;++i)    
      if (CPU_ISSET(i, &cpuset))
	  std::cout << i << std::endl;
    */
    
#endif
}

void bench(unsigned int threads, const std::string& func_name, double dx, double dy, void (*bench_func)(double ddx, double ddy, long long& ret))
{
    std::cout << "Start " << func_name << " " << threads << std::endl;
    auto start = get_time(0);

    auto v_threads = std::vector<std::thread>(threads);
    auto totals = std::vector<long long>(threads);
    for (unsigned int i = 0; i < threads; i++)
    {
        v_threads[i] = std::thread(bench_func, dx, dy, std::ref(totals[i]));
        if (G_pin_threads)
        {
            const auto thread_handle = v_threads[i].native_handle();
            const unsigned int ht_thread_per_core = 2;
            const unsigned int max_cores_non_ht = std::thread::hardware_concurrency() / ht_thread_per_core;
            const unsigned int core_id = i < max_cores_non_ht ? ht_thread_per_core * i : ht_thread_per_core * (i - max_cores_non_ht) + 1; // first real cores, then HT
            set_thread_affinity(thread_handle, core_id);
        }
    }
    for (unsigned int i = 0; i < threads; i++)
        v_threads[i].join();
    long long total = 0;
    for (unsigned int i = 0; i < threads; i++)
        total += totals[i];
    auto end = get_time(total);
    std::cout << "res=" << total << " " << total / threads << " " << total % threads << std::endl;
    std::chrono::duration<double> diff = end - start;
    std::cout << diff.count() << " s" << " " << diff.count() / threads << " s/thread" << std::endl;
    std::cout << "==============" << std::endl;
}

void bench_threads(const std::string& func_name, double dx, double dy, void (*bench_func)(double ddx, double ddy, long long& ret))
{
    for (unsigned int threads = 1; threads <= std::thread::hardware_concurrency(); threads++)
        bench(threads, func_name, dx, dy, bench_func);
}

void test_polynom()
{
    double res;
    double maximum = -1;
    for (int s = 0; s < 5000; s++)
    {
        std::cout << s;
        std::cout.flush();
        polynom2 p(s);
        res = p.find_roots(150);
        std::cout << " : " << res << std::endl;
        maximum = std::max(maximum, res);
    }
    std::cout << "max=" << maximum << std::endl;
}

void test_polynom0()
{
    int s = 656;
    polynom2 p(s);
    auto res = p.find_roots(150);
    std::cout << "res" << s << ": " << res << std::endl << std::endl;
}

// Warning: feature check will not run if specific instructions are used for prologue
int effective_main(int argc, char** argv)
{
    std::cout << "Starting..." << std::endl;

    //test_polynom(); return 0;

#if defined(_MSC_FULL_VER)
    std::cout << "MSC_FULL_VER: " << _MSC_FULL_VER << std::endl;
#endif

    // Cf https://docs.microsoft.com/fr-fr/cpp/preprocessor/predefined-macros?view=msvc-160
#if defined(__AVX__)
    std::cout << "__AVX__: " << __AVX__ << std::endl;
    if (!check_avx())
        std::cerr << "Avx seems not supported!" << std::endl;
#endif

#if defined(__AVX2__)
    std::cout << "__AVX2__: " << __AVX2__ << std::endl;
    if (!check_avx2())
        std::cerr << "Avx2 seems not supported!" << std::endl;
#endif

#if defined(__AVX512F__)
    std::cout << "__AVX512F__: " << __AVX512F__ << std::endl;
    if (!check_avx512())
        std::cerr << "Avx512 seems not supported!" << std::endl;
#endif

#if defined(__GNUC__) 
    std::cout << "g++ compiler: " << __GNUC__ << "." << __GNUC_MINOR__ << std::endl;
#endif

#if defined (__clang_major__)
    std::cout << "Clang compiler: " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__ << std::endl;
#endif

#if defined (__INTEL_LLVM_COMPILER)
    std::cout << "Intel LLVM compiler: " << __VERSION__ << std::endl;
#endif	

    if (argc == 1)
    {
        std::cout << "Syntax: " << argv[0] << " (Mandel, Mandel2, Heat, Pi, Fft, Matrix, Polynom, All)" << std::endl;
        return 0;
    }

    bool doAll = (strcmp(argv[1], "All") == 0);
    for (int i = 1; i < argc; i++)
    {
        if ((strcmp(argv[i], "Mandel") == 0) || doAll)
            bench_threads("benchMandel", 0.0005, 0.0005, benchMandel);
        if ((strcmp(argv[i], "Mandel2") == 0) || doAll)
            bench_threads("benchMandel2", 0.0005, 0.0005, benchMandel2);
        if ((strcmp(argv[i], "Heat") == 0) || doAll)
            bench_threads("benchHeat", 0.01, 0.01, benchHeat);
        if ((strcmp(argv[i], "Pi") == 0) || doAll)
            bench_threads("benchPi", 1e-9, 0.0, benchPi);
        if ((strcmp(argv[i], "Fft") == 0) || doAll)
            bench_threads("benchFft", 24, 0.1234, benchFft);
        if ((strcmp(argv[i], "Matrix") == 0) || doAll)
            bench_threads("benchMatrix", 1000, 0.0, benchMatrix);
        if ((strcmp(argv[i], "Polynom") == 0) || doAll)
            bench_threads("benchPolynom", 150, 5000, benchPolynom);
    }
    return 0;
}

void parseOption(int argc, char** argv, bool& throwFPE, bool& pinThreads)
{
    for (int i = 0; i < argc; i++)
    {
        if (argv[i] == std::string_view("-throwfpe"))
            throwFPE = true;
        if (argv[i] == std::string_view("-pinthreads"))
            pinThreads = true;
    }
}

#if defined(_WIN32)
void setThrowFPE()
{
    unsigned int cw;
    _controlfp_s(&cw, 0, 0);
    unsigned int new_value = cw & ~(_EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW);
    _controlfp_s(&cw, new_value, _MCW_EM);
    std::cout << "**** Throw FPE activated ****" << std::endl;
}

int filter_exception(unsigned int code, struct _EXCEPTION_POINTERS* ep)
{
    std::cout << std::format("Exception code: 0x{:x}", code);
    if (ep != nullptr && ep->ExceptionRecord != nullptr)
        std::cout << ", Address=0x" << (void*)ep->ExceptionRecord->ExceptionAddress;
    std::cout << " : ";
    switch (code)
    {
    case  EXCEPTION_ACCESS_VIOLATION:
        std::cout << "Acess Violation" << std::endl;
        break;
    case  EXCEPTION_DATATYPE_MISALIGNMENT:
        std::cout << "Data Type MisAlignement" << std::endl;
        break;
    case  EXCEPTION_BREAKPOINT:
        std::cout << "Breakpoint" << std::endl;
        break;
    case  EXCEPTION_SINGLE_STEP:
        std::cout << "Single Step" << std::endl;
        break;
    case  EXCEPTION_ARRAY_BOUNDS_EXCEEDED:
        std::cout << "Array Bounds Exceeded" << std::endl;
        break;
    case  EXCEPTION_FLT_DENORMAL_OPERAND:
        std::cout << "Floating Point Denormal Operand" << std::endl;
        break;
    case  EXCEPTION_FLT_DIVIDE_BY_ZERO:
        std::cout << "Floating Point Divide By Zero" << std::endl;
        break;
    case  EXCEPTION_FLT_INEXACT_RESULT:
        std::cout << "Floating Point Inexact Result" << std::endl;
        break;
    case  EXCEPTION_FLT_INVALID_OPERATION:
        std::cout << "Floating Point Invalid Operation" << std::endl;
        break;
    case  EXCEPTION_FLT_OVERFLOW:
        std::cout << "Floating Point Overflow" << std::endl;
        break;
    case  EXCEPTION_FLT_STACK_CHECK:
        std::cout << "Floating Point Stack Check" << std::endl;
        break;
    case  EXCEPTION_FLT_UNDERFLOW:
        std::cout << "Floating Point Underflow" << std::endl;
        break;
    case  EXCEPTION_INT_DIVIDE_BY_ZERO:
        std::cout << "Integer Divide By Zero" << std::endl;
        break;
    case  EXCEPTION_INT_OVERFLOW:
        std::cout << "Integer Overflow" << std::endl;
        break;
    case  EXCEPTION_PRIV_INSTRUCTION:
        std::cout << "Privileged Instruction" << std::endl;
        break;
    case  EXCEPTION_IN_PAGE_ERROR:
        std::cout << "Page Error" << std::endl;
        break;
    case  EXCEPTION_ILLEGAL_INSTRUCTION:
        std::cout << "Illegal Instruction" << std::endl;
        break;
    case  EXCEPTION_NONCONTINUABLE_EXCEPTION:
        std::cout << "Non Continuable Exception" << std::endl;
        break;
    case  EXCEPTION_STACK_OVERFLOW:
        std::cout << "Stack Overflow" << std::endl;
        break;
    case  EXCEPTION_INVALID_DISPOSITION:
        std::cout << "Invalid Disposition" << std::endl;
        break;
    case  EXCEPTION_GUARD_PAGE:
        std::cout << "Guard Page" << std::endl;
        break;
    case  EXCEPTION_INVALID_HANDLE:
        std::cout << "Invalid Handle" << std::endl;
        break;
    case  CONTROL_C_EXIT:
        std::cout << "Control-C Exit" << std::endl;
        break;
    default:
        std::cout << "Unknown Exception" << std::endl;
        break;
    }
    return EXCEPTION_EXECUTE_HANDLER;
}
#endif

int main(int argc, char** argv)
{
    bool throwFPE = false;
    G_pin_threads = false;
    parseOption(argc, argv, throwFPE, G_pin_threads);
    if (G_pin_threads)
        std::cout << "** Pin threads activated" << std::endl;
#if defined(_WIN32)
    __try
    {
        if (throwFPE)
            setThrowFPE();
        /*
               double a = 1.0;
               double b = 0.0;
               double c = a/b;
               std::cout << c << std::endl;
        */

        effective_main(argc, argv);
    }
    __except (filter_exception(GetExceptionCode(), GetExceptionInformation()))
    {
        exit(GetExceptionCode());
    }
#else
    effective_main(argc, argv);
#endif
}
