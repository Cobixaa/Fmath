#ifndef FMATH_H
#define FMATH_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Configuration
#ifndef FMATH_TABLE_BITS
#define FMATH_TABLE_BITS 12 /* 4096-entry LUT by default */
#endif

#ifndef FMATH_LUT_INIT_WITH_LIBM
#define FMATH_LUT_INIT_WITH_LIBM 1 /* use sinf/cosf from libm at init */
#endif

#ifndef FMATH_ENABLE_OMP
#define FMATH_ENABLE_OMP 0
#endif

#if FMATH_ENABLE_OMP
#include <omp.h>
#endif

#ifndef FMATH_ALWAYS_INLINE
#if defined(__GNUC__) || defined(__clang__)
#define FMATH_ALWAYS_INLINE __attribute__((always_inline))
#else
#define FMATH_ALWAYS_INLINE
#endif
#endif

#ifndef FMATH_INLINE
#define FMATH_INLINE static inline FMATH_ALWAYS_INLINE
#endif

// Public API

// Initializes lookup tables (called lazily by first use). Safe to call multiple times.
void fmath_init(void);

// Scalar fast approximations (single-precision)
float fmath_sinf(float x);
float fmath_cosf(float x);
float fmath_expf(float x);
float fmath_logf(float x);
float fmath_sqrtf(float x);
float fmath_rsqrtf(float x);
float fmath_rcpf(float x);

// Array APIs (in-place allowed if dst == src)
void fmath_sinf_array(float *dst, const float *src, size_t count);
void fmath_cosf_array(float *dst, const float *src, size_t count);
void fmath_expf_array(float *dst, const float *src, size_t count);
void fmath_logf_array(float *dst, const float *src, size_t count);
void fmath_sqrtf_array(float *dst, const float *src, size_t count);
void fmath_rsqrtf_array(float *dst, const float *src, size_t count);
void fmath_rcpf_array(float *dst, const float *src, size_t count);

#ifdef __cplusplus
}
#endif

#endif /* FMATH_H */