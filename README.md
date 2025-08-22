fmath: Fast Approximate Math for C (float32)

Overview
--------
`fmath` is a high-performance C library providing ultra-fast approximations for common math functions: sin, cos, exp, log, sqrt, rsqrt, and reciprocal. It prioritizes throughput over strict IEEE accuracy, targeting single-precision floats by default.

No Makefile — Pure Command Build
--------------------------------
Build and run the benchmark with a single command (Linux, gcc/clang):

```bash
gcc -O3 -ffast-math -march=native -funroll-loops -Iinclude src/fmath.c bench/bench.c -o fmath_bench -lm
./fmath_bench            # optional: ./fmath_bench <N>, default N=8000000
```

Enable OpenMP (optional):

```bash
gcc -O3 -ffast-math -march=native -funroll-loops -fopenmp -DFMATH_ENABLE_OMP=1 -Iinclude src/fmath.c bench/bench.c -o fmath_bench -lm
```

Easy API
--------
Short names: define `FMATH_SHORT_NAMES` before including `fmath.h`.

```c
#define FMATH_SHORT_NAMES 1
#include "fmath.h"

void demo() {
	fm_init();
	float x = 1.2f;
	float s = fm_sin(x);
	float c = fm_cos(x);
	float e = fm_exp(x);
	float l = fm_log(x);
	float r = fm_rsqrt(3.0f);
}
```

Override libm (opt-in):

```c
#define FMATH_OVERRIDE_LIBM 1
#include "fmath.h"   // sinf/cosf/expf/logf/sqrtf map to fmath versions
```

Array helpers:
- `fmath_*_array(dst, src, count)` process arrays
- With `FMATH_SHORT_NAMES`: `fm_*_arr(dst, src, n)` and `fm_*_aa(dst, src)` (count-of src)

What’s Implemented (Fast Paths)
-------------------------------
- `sin, cos`: LUT (2^FMATH_TABLE_BITS, default 4096) + linear interpolation; `cos` via phase shift
- `exp`: magic-bias range reduction r=x*log2(e)=n+f; cubic for 2^f; scale by 2^n via exponent bits
- `log`: extract exponent/mantissa; 5-term `log(1+z)` polynomial
- `rsqrt`: Quake constant + 1 Newton step
- `sqrt`: `x * rsqrt(x)`
- `rcp`: `1/x` (can be swapped for NR refine if desired)

Tuning and Options
------------------
- `FMATH_TABLE_BITS` (default 12): LUT size for sin/cos
- `FMATH_LUT_INIT_WITH_LIBM` (default 1): init LUT via `sinf`
- `FMATH_ENABLE_OMP` (default 0): enable OpenMP for array functions
- `FMATH_SHORT_NAMES`: short API aliases
- `FMATH_OVERRIDE_LIBM`: redefine `sinf/cosf/expf/logf/sqrtf` to fmath variants

Faster Than libm — Notes
------------------------
- The provided benchmark shows substantial speedups for sin/cos/log/sqrt/rsqrt/rcp.
- `exp` now uses a very fast range-reduction with a short cubic and is competitive or faster. If your CPU/libc combination still outperforms, try:
  - Add `-mfma` if supported; or build with `-ffp-contract=fast` (implicit with `-ffast-math`).
  - Keep inputs within typical ranges (e.g., [-10, 10]) for best accuracy and speed.

Tutorial: Fast Math Techniques
------------------------------
1) Lookup Tables (LUTs)
- Precompute sin/cos samples, use linear or cubic interpolation
- O(1) speed, higher memory usage

2) Polynomial Approximations
- Taylor, Chebyshev, minimax; small-degree polynomials for core ranges

3) Rational Approximations (Padé)
- Ratio of polynomials; often more accurate than plain polynomials

4) Range Reduction
- Bring inputs into a small interval where approximations are accurate

5) Bit Hacks
- Manipulate IEEE-754 bits (e.g., rsqrt, log/exp exponent/mantissa)

6) Piecewise Approximations
- Split domain into intervals with dedicated polynomials per segment

7) CORDIC Algorithm
- Shift/add iterative trig/log without multiplies (great for embedded)

8) Horner’s Rule
- Evaluate polynomials efficiently to reduce multiplies

9) Vectorization (SIMD)
- SSE/AVX/NEON to process many floats in parallel

10) Loop Unrolling & Pipelining
- Reduce branch overhead, keep execution units busy

11) Fused Multiply-Add (FMA)
- `a*b + c` in one instruction; faster and more stable

12) Avoid Denormals
- Add tiny constant or set FTZ/DAZ to avoid slow paths on denormals

13) Pre-Normalization
- Normalize/scaling inputs to keep approximations stable

14) Reduced Precision
- Prefer float32 when acceptable for 2x speed and half bandwidth

15) Speculative Execution
- Compute candidates and select mask-wise to avoid branches

16) Table Compression
- Store LUT in smaller integer formats and rescale

17) Hybrid Methods
- Mix LUT and polynomials for optimal memory vs speed

18) Recurrence Relations
- Use math identities to compute sequences cheaply

19) Inlining
- `static inline` hot paths to remove call overhead

20) Compiler Flags
- `-O3 -ffast-math -march=native -funroll-loops` for aggressive optimization

21) Parallelization
- Use threads/OpenMP/GPU for array workloads

22) Approximate Reciprocal & Division
- Hardware reciprocal + Newton step to refine

23) Branchless Tricks
- Replace conditionals with masks to avoid mispredictions

24) Pre-Transform Inputs
- Pre-scale or pre-index to save conversions in hot loops

25) Function Specialization
- Handle frequent special cases early (e.g., `exp(0)=1`)

Benchmarks
----------
Run the included benchmark to compare `fmath` vs `libm` on large arrays. Results vary by CPU, compiler, and flags, but the aim is for every function to be faster than stock libm in hot loops.

Accuracy Notes
--------------
- Approximations trade precision for speed; test against your workload and add refinements (e.g., second Newton step for `rsqrt`) when needed.

License
-------
MIT