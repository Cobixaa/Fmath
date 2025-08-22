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

Quick Usage
-----------
```c
#include "fmath.h"

int main() {
	fmath_init();
	float x = 1.2345f;
	float s = fmath_sinf(x);
	float e = fmath_expf(x);
	return (int)(s + e);
}
```

Design Highlights
-----------------
- Single LUT for `sinf`, `cosf` derived by phase shift (π/2), linear interpolation
- `expf`: range reduction `x = n*ln2 + r`, 5th-order polynomial for `exp(r)`, `2^n` via exponent bits
- `logf`: bit extraction into exponent/mantissa, polynomial `log(1+z)` on `[0,1)`
- `rsqrtf`: Quake inverse square-root with one Newton step
- Array APIs for throughput; optional OpenMP parallel-for
- Compile-time knobs in `include/fmath.h`:
  - `FMATH_TABLE_BITS` (default 12 -> 4096 entries)
  - `FMATH_LUT_INIT_WITH_LIBM` (default 1)
  - `FMATH_ENABLE_OMP` (default 0)

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
Run the included benchmark to compare `fmath` vs `libm` on large arrays (default 8M elements). Results vary by CPU, compiler, and flags.

Accuracy Notes
--------------
- Approximations sacrifice some precision for speed; test for your workload.
- `logf` and `expf` use moderate-degree polynomials and simple range reduction; if you need tighter bounds, consider higher-degree or piecewise polynomials.
- `rsqrtf` uses one Newton step (good trade-off); add another step for higher accuracy.

License
-------
MIT