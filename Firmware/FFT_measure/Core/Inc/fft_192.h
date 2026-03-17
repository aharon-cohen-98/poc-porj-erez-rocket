/**
  ******************************************************************************
  * @file    fft_192.h
  * @brief   192-point FFT via 3x64 composite using CMSIS-DSP 64-point CFFT
  ******************************************************************************
  */

#ifndef FFT_192_H
#define FFT_192_H

#ifdef __cplusplus
extern "C" {
#endif

#include "arm_math.h"

/** Number of real samples for 192-point real FFT input */
#define FFT_192_INPUT_LEN  192U

/** Number of floats for 192-point complex output (192 * 2 = 384) */
#define FFT_192_OUTPUT_LEN 384U

/**
  * @brief  192-point FFT from real input using 3x64 composite (64-point CFFT).
  * @param  x_192   Input: 192 real samples
  * @param  y_384   Output: 192 complex bins, interleaved [Re0,Im0, Re1,Im1, ...] (384 floats)
  * @param  cfft_64  Initialized 64-point complex FFT instance (arm_cfft_init_f32(S, 64))
  * @note   For real input, only bins 0..96 are unique (conjugate symmetry: bin k = conj(bin 192-k)).
  *         Bins 0 and 96 are real; use y_384[0], y_384[192] for DC and Nyquist.
  */
void fft_192_real_f32(const float32_t *x_192, float32_t *y_384,
                      const arm_cfft_instance_f32 *cfft_64);

/**
  * @brief  192-point FFT from complex input using 3x64 composite (64-point CFFT).
  * @param  x_384   Input: 192 complex samples, interleaved [Re0,Im0, Re1,Im1, ...] (384 floats)
  * @param  y_384   Output: 192 complex bins, interleaved (384 floats)
  * @param  cfft_64  Initialized 64-point complex FFT instance
  */
void fft_192_cplx_f32(const float32_t *x_384, float32_t *y_384,
                      const arm_cfft_instance_f32 *cfft_64);

#ifdef __cplusplus
}
#endif

#endif /* FFT_192_H */
