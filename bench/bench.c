#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#include "fmath.h"

static double now_time(void) {
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static uint32_t rng_state = 0x12345678u;
static inline uint32_t xorshift32(void) {
	uint32_t x = rng_state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return rng_state = x;
}

static inline float randf_range(float a, float b) {
	uint32_t r = xorshift32();
	float t = (float)(r >> 8) * (1.0f / 16777216.0f); // 24-bit mantissa
	return a + (b - a) * t;
}

static void fill_range(float *arr, size_t n, float a, float b) {
	for (size_t i = 0; i < n; ++i) arr[i] = randf_range(a, b);
}

static void fill_positive(float *arr, size_t n, float minv, float maxv) {
	for (size_t i = 0; i < n; ++i) {
		float v = randf_range(minv, maxv);
		if (v <= 0.0f) v = minv;
		arr[i] = v;
	}
}

static void fill_nonzero(float *arr, size_t n, float min_abs, float max_abs) {
	for (size_t i = 0; i < n; ++i) {
		float v = randf_range(-max_abs, max_abs);
		if (v >= 0.0f && v < min_abs) v = min_abs;
		if (v < 0.0f && v > -min_abs) v = -min_abs;
		arr[i] = v;
	}
}

static double time_loop(float *dst, const float *src, size_t n, float (*fn)(float)) {
	double t0 = now_time();
	for (size_t i = 0; i < n; ++i) dst[i] = fn(src[i]);
	return now_time() - t0;
}

static float rsqrt_libm(float x) {
	if (x <= 0.0f) return NAN;
	return 1.0f / sqrtf(x);
}

static float rcp_libm(float x) {
	if (x == 0.0f) return copysignf(INFINITY, x);
	return 1.0f / x;
}

int main(int argc, char **argv) {
	size_t n = (argc > 1) ? (size_t)atoll(argv[1]) : (size_t)8 * 1000 * 1000; // default 8M
	printf("fmath bench n=%zu\n", n);

	float *in = (float*)malloc(n * sizeof(float));
	float *out = (float*)malloc(n * sizeof(float));
	if (!in || !out) {
		fprintf(stderr, "allocation failed\n");
		return 1;
	}

	fmath_init();

	// sin
	fill_range(in, n, -1000.0f, 1000.0f);
	double t_fmath = now_time();
	fmath_sinf_array(out, in, n);
	t_fmath = now_time() - t_fmath;
	double t_libm = time_loop(out, in, n, sinf);
	printf("sin: fmath=%.3f s, libm=%.3f s, speedup=%.2fx\n", t_fmath, t_libm, t_libm / t_fmath);

	// cos
	fill_range(in, n, -1000.0f, 1000.0f);
	t_fmath = now_time();
	fmath_cosf_array(out, in, n);
	t_fmath = now_time() - t_fmath;
	t_libm = time_loop(out, in, n, cosf);
	printf("cos: fmath=%.3f s, libm=%.3f s, speedup=%.2fx\n", t_fmath, t_libm, t_libm / t_fmath);

	// exp
	fill_range(in, n, -10.0f, 10.0f);
	t_fmath = now_time();
	fmath_expf_array(out, in, n);
	t_fmath = now_time() - t_fmath;
	t_libm = time_loop(out, in, n, expf);
	printf("exp: fmath=%.3f s, libm=%.3f s, speedup=%.2fx\n", t_fmath, t_libm, t_libm / t_fmath);

	// log
	fill_positive(in, n, 1e-6f, 1e6f);
	t_fmath = now_time();
	fmath_logf_array(out, in, n);
	t_fmath = now_time() - t_fmath;
	t_libm = time_loop(out, in, n, logf);
	printf("log: fmath=%.3f s, libm=%.3f s, speedup=%.2fx\n", t_fmath, t_libm, t_libm / t_fmath);

	// sqrt
	fill_positive(in, n, 1e-6f, 1e6f);
	t_fmath = now_time();
	fmath_sqrtf_array(out, in, n);
	t_fmath = now_time() - t_fmath;
	t_libm = time_loop(out, in, n, sqrtf);
	printf("sqrt: fmath=%.3f s, libm=%.3f s, speedup=%.2fx\n", t_fmath, t_libm, t_libm / t_fmath);

	// rsqrt
	fill_positive(in, n, 1e-6f, 1e6f);
	t_fmath = now_time();
	fmath_rsqrtf_array(out, in, n);
	t_fmath = now_time() - t_fmath;
	t_libm = time_loop(out, in, n, rsqrt_libm);
	printf("rsqrt: fmath=%.3f s, libm=%.3f s, speedup=%.2fx\n", t_fmath, t_libm, t_libm / t_fmath);

	// rcp
	fill_nonzero(in, n, 1e-3f, 1e6f);
	t_fmath = now_time();
	fmath_rcpf_array(out, in, n);
	t_fmath = now_time() - t_fmath;
	t_libm = time_loop(out, in, n, rcp_libm);
	printf("rcp: fmath=%.3f s, libm=%.3f s, speedup=%.2fx\n", t_fmath, t_libm, t_libm / t_fmath);

	free(in);
	free(out);
	return 0;
}