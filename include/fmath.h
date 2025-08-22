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

// Easy API helpers
#ifndef FMATH_COUNT_OF
#define FMATH_COUNT_OF(arr) (sizeof(arr) / sizeof((arr)[0]))
#endif

#ifdef FMATH_SHORT_NAMES
#define fm_init            fmath_init
#define fm_sin             fmath_sinf
#define fm_cos             fmath_cosf
#define fm_exp             fmath_expf
#define fm_log             fmath_logf
#define fm_sqrt            fmath_sqrtf
#define fm_rsqrt           fmath_rsqrtf
#define fm_rcp             fmath_rcpf
#define fm_sin_arr(dst, src, n)   fmath_sinf_array((dst), (src), (n))
#define fm_cos_arr(dst, src, n)   fmath_cosf_array((dst), (src), (n))
#define fm_exp_arr(dst, src, n)   fmath_expf_array((dst), (src), (n))
#define fm_log_arr(dst, src, n)   fmath_logf_array((dst), (src), (n))
#define fm_sqrt_arr(dst, src, n)  fmath_sqrtf_array((dst), (src), (n))
#define fm_rsqrt_arr(dst, src, n) fmath_rsqrtf_array((dst), (src), (n))
#define fm_rcp_arr(dst, src, n)   fmath_rcpf_array((dst), (src), (n))
#define fm_sin_aa(dst, src)       fmath_sinf_array((dst), (src), FMATH_COUNT_OF(src))
#define fm_cos_aa(dst, src)       fmath_cosf_array((dst), (src), FMATH_COUNT_OF(src))
#define fm_exp_aa(dst, src)       fmath_expf_array((dst), (src), FMATH_COUNT_OF(src))
#define fm_log_aa(dst, src)       fmath_logf_array((dst), (src), FMATH_COUNT_OF(src))
#define fm_sqrt_aa(dst, src)      fmath_sqrtf_array((dst), (src), FMATH_COUNT_OF(src))
#define fm_rsqrt_aa(dst, src)     fmath_rsqrtf_array((dst), (src), FMATH_COUNT_OF(src))
#define fm_rcp_aa(dst, src)       fmath_rcpf_array((dst), (src), FMATH_COUNT_OF(src))
#endif

#ifdef FMATH_OVERRIDE_LIBM
// Override common libm float entry points with fmath versions
// Toggle by defining FMATH_OVERRIDE_LIBM before including this header
#ifndef FMATH_NO_OVERRIDE_SIN
#define sinf fmath_sinf
#endif
#ifndef FMATH_NO_OVERRIDE_COS
#define cosf fmath_cosf
#endif
#ifndef FMATH_NO_OVERRIDE_EXP
#define expf fmath_expf
#endif
#ifndef FMATH_NO_OVERRIDE_LOG
#define logf fmath_logf
#endif
#ifndef FMATH_NO_OVERRIDE_SQRT
#define sqrtf fmath_sqrtf
#endif
#endif

#ifdef __cplusplus
}
#endif

#endif /* FMATH_H */