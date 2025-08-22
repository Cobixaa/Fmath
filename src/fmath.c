#include "fmath.h"

#include <math.h>
#include <string.h>
#include <stdbool.h>

#ifndef FMATH_PI
#define FMATH_PI 3.14159265358979323846f
#endif
#ifndef FMATH_TWO_PI
#define FMATH_TWO_PI 6.28318530717958647692f
#endif
#ifndef FMATH_INV_TWO_PI
#define FMATH_INV_TWO_PI 0.15915494309189533577f /* 1/(2*pi) */
#endif
#ifndef FMATH_LN2
#define FMATH_LN2 0.69314718055994530942f
#endif
#ifndef FMATH_INV_LN2
#define FMATH_INV_LN2 1.4426950408889634074f /* 1/ln(2) */
#endif

// Internal LUT for sin; cos derived via phase shift
enum {
	FMATH_TABLE_SIZE = 1 << FMATH_TABLE_BITS,
	FMATH_TABLE_MASK = FMATH_TABLE_SIZE - 1
};

static float fmath_sin_lut[FMATH_TABLE_SIZE];
static float fmath_index_scale = 0.0f; // FMATH_TABLE_SIZE / (2*pi)
static int fmath_is_initialized = 0;

FMATH_INLINE void fmath_init_once(void) {
	if (fmath_is_initialized) return;
	fmath_index_scale = (float)FMATH_TABLE_SIZE * (1.0f / FMATH_TWO_PI);
	const float step = FMATH_TWO_PI / (float)FMATH_TABLE_SIZE;
#if FMATH_LUT_INIT_WITH_LIBM
	for (int i = 0; i < FMATH_TABLE_SIZE; ++i) {
		fmath_sin_lut[i] = sinf(step * (float)i);
	}
#else
	// Fallback: compute via 5th-order Taylor around 0 (less accurate)
	for (int i = 0; i < FMATH_TABLE_SIZE; ++i) {
		float x = step * (float)i;
		float x2 = x * x;
		float x3 = x2 * x;
		float x5 = x3 * x2;
		fmath_sin_lut[i] = x - (x3 * (1.0f / 6.0f)) + (x5 * (1.0f / 120.0f));
	}
#endif
	fmath_is_initialized = 1;
}

void fmath_init(void) {
	fmath_init_once();
}

// Helpers for bit-casting without aliasing UB
FMATH_INLINE uint32_t fmath_bitcast_f32_to_u32(float x) {
	uint32_t u;
	memcpy(&u, &x, sizeof u);
	return u;
}

FMATH_INLINE float fmath_bitcast_u32_to_f32(uint32_t u) {
	float x;
	memcpy(&x, &u, sizeof x);
	return x;
}

// Fast sinf/cosf using LUT + linear interpolation, with power-of-two table size.
float fmath_sinf(float x) {
	if (!fmath_is_initialized) fmath_init_once();
	// Map x radians to table index space, wrapping via mask
	float index_f = x * fmath_index_scale;
	float idx_floor = floorf(index_f);
	int i0 = ((int)idx_floor) & FMATH_TABLE_MASK;
	int i1 = (i0 + 1) & FMATH_TABLE_MASK;
	float t = index_f - idx_floor;
	float s0 = fmath_sin_lut[i0];
	float s1 = fmath_sin_lut[i1];
	return s0 + t * (s1 - s0);
}

float fmath_cosf(float x) {
	if (!fmath_is_initialized) fmath_init_once();
	// cos(x) = sin(x + pi/2) -> phase shift by quarter table
	float index_f = (x + 0.5f * FMATH_PI) * fmath_index_scale;
	float idx_floor = floorf(index_f);
	int i0 = ((int)idx_floor) & FMATH_TABLE_MASK;
	int i1 = (i0 + 1) & FMATH_TABLE_MASK;
	float t = index_f - idx_floor;
	float s0 = fmath_sin_lut[i0];
	float s1 = fmath_sin_lut[i1];
	return s0 + t * (s1 - s0);
}

// Fast expf using range reduction x = n*ln2 + r, |r| <= ln2/2, and 5th-order poly
float fmath_expf(float x) {
	// Clamp to avoid extreme overflow behavior
	if (x > 88.0f) return INFINITY;        // ~expf overflow threshold for float
	if (x < -100.0f) return 0.0f;          // underflow to 0

	int n = (int)(x * FMATH_INV_LN2 + (x >= 0.0f ? 0.5f : -0.5f));
	float r = x - (float)n * FMATH_LN2;
	// Polynomial: exp(r) ≈ 1 + r + r^2/2 + r^3/6 + r^4/24 + r^5/120
	float p = 1.0f + r * (1.0f + r * (0.5f + r * (0.16666667163372f + r * (0.04166666790843f + r * 0.00833333376795f))));
	// Compute 2^n quickly when in normal range; fallback to ldexpf otherwise
	float two_n;
	if (n >= -126 && n <= 127) {
		uint32_t bits = (uint32_t)(n + 127) << 23;
		two_n = fmath_bitcast_u32_to_f32(bits);
	} else {
		two_n = ldexpf(1.0f, n);
	}
	return p * two_n;
}

// Fast logf using bit tricks: x = m * 2^e with m in [1,2). log(x)=e*ln2 + log(m)
float fmath_logf(float x) {
	if (x <= 0.0f) {
		if (x == 0.0f) return -INFINITY;
		return NAN;
	}
	uint32_t xi = fmath_bitcast_f32_to_u32(x);
	int e = (int)((xi >> 23) & 255) - 127;
	uint32_t mant = (xi & 0x7fffffU) | 0x3f800000U; // 1.m
	float m = fmath_bitcast_u32_to_f32(mant);
	float z = m - 1.0f; // in [0,1)
	// log(1+z) via alternating series with 5 terms: z - z^2/2 + z^3/3 - z^4/4 + z^5/5
	float z2 = z * z;
	float z3 = z2 * z;
	float z4 = z3 * z;
	float z5 = z4 * z;
	float log1pz = z - 0.5f * z2 + (z3 * 0.3333333433f) - 0.25f * z4 + 0.2f * z5;
	return (float)e * FMATH_LN2 + log1pz;
}

// Fast inverse sqrt (Quake III) + one Newton-Raphson refinement
float fmath_rsqrtf(float x) {
	if (x <= 0.0f) {
		if (x == 0.0f) return INFINITY;
		return NAN;
	}
	float xhalf = 0.5f * x;
	uint32_t i = fmath_bitcast_f32_to_u32(x);
	i = 0x5f3759dfu - (i >> 1);
	float y = fmath_bitcast_u32_to_f32(i);
	// One NR iteration
	y = y * (1.5f - xhalf * y * y);
	return y;
}

float fmath_sqrtf(float x) {
	if (x <= 0.0f) {
		if (x == 0.0f) return 0.0f;
		return NAN;
	}
	return x * fmath_rsqrtf(x);
}

float fmath_rcpf(float x) {
	if (x == 0.0f) {
#if defined(__GNUC__) || defined(__clang__)
		return copysignf(INFINITY, x);
#else
		return (x < 0.0f) ? -INFINITY : INFINITY;
#endif
	}
	// Fallback to hardware division; can be replaced with NR refine if desired
	return 1.0f / x;
}

// Array APIs
void fmath_sinf_array(float *dst, const float *src, size_t count) {
	if (!fmath_is_initialized) fmath_init_once();
	#if FMATH_ENABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (size_t i = 0; i < count; ++i) dst[i] = fmath_sinf(src[i]);
}

void fmath_cosf_array(float *dst, const float *src, size_t count) {
	if (!fmath_is_initialized) fmath_init_once();
	#if FMATH_ENABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (size_t i = 0; i < count; ++i) dst[i] = fmath_cosf(src[i]);
}

void fmath_expf_array(float *dst, const float *src, size_t count) {
	#if FMATH_ENABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (size_t i = 0; i < count; ++i) dst[i] = fmath_expf(src[i]);
}

void fmath_logf_array(float *dst, const float *src, size_t count) {
	#if FMATH_ENABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (size_t i = 0; i < count; ++i) dst[i] = fmath_logf(src[i]);
}

void fmath_sqrtf_array(float *dst, const float *src, size_t count) {
	#if FMATH_ENABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (size_t i = 0; i < count; ++i) dst[i] = fmath_sqrtf(src[i]);
}

void fmath_rsqrtf_array(float *dst, const float *src, size_t count) {
	#if FMATH_ENABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (size_t i = 0; i < count; ++i) dst[i] = fmath_rsqrtf(src[i]);
}

void fmath_rcpf_array(float *dst, const float *src, size_t count) {
	#if FMATH_ENABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (size_t i = 0; i < count; ++i) dst[i] = fmath_rcpf(src[i]);
}