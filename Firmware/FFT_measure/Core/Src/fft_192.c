/**
  ******************************************************************************
  * @file    fft_192.c
  * @brief   192-point FFT via 3x64 composite (Cooley-Tukey): 64x 3-point DFT,
  *          twiddle multiply, then 3x 64-point arm_cfft_f32.
  ******************************************************************************
  */

#include "fft_192.h"
#include <stdint.h>
#include <string.h>

#define N1       3U
#define N2      64U
#define N_FFT  192U

/* Working buffer: 3 rows x 64 complex = 384 floats, row-major (each row 128 floats) */
__attribute__((aligned(32))) static float32_t work_3x64[N1 * N2 * 2];

/* Twiddle W_192^(64*k1*n2) = W3^(k1*n2). Only 3 unique values: 1, W3, W3^2 */
static float32_t tw3_re[3], tw3_im[3];
static uint8_t tw3_init;

static void init_twiddles(void)
{
  if (tw3_init)
    return;
  /* W3 = exp(-j*2*pi/3) = cos(2*pi/3) - j*sin(2*pi/3) = -0.5 - j*sqrt(3)/2 */
  const float32_t two_pi = 2.0f * 3.14159265358979323846f;
  tw3_re[0] = 1.0f;
  tw3_im[0] = 0.0f;
  tw3_re[1] = -0.5f;
  tw3_im[1] = -0.86602540378443864676f;  /* -sqrt(3)/2 */
  tw3_re[2] = -0.5f;
  tw3_im[2] =  0.86602540378443864676f;  /* +sqrt(3)/2 */
  tw3_init = 1;
}

/* 3-point complex DFT: Y(k) = sum_{n=0..2} x(n)*W3^(k*n). In/out interleaved [Re,Im,...] */
static void dft3_cplx_f32(const float32_t *x, float32_t *y)
{
  const float32_t x0r = x[0], x0i = x[1];
  const float32_t x1r = x[2], x1i = x[3];
  const float32_t x2r = x[4], x2i = x[5];

  /* Y(0) = x0 + x1 + x2 */
  y[0] = x0r + x1r + x2r;
  y[1] = x0i + x1i + x2i;

  /* Y(1) = x0 + W3*x1 + W3^2*x2 */
  {
    const float32_t w1r = tw3_re[1], w1i = tw3_im[1];
    const float32_t w2r = tw3_re[2], w2i = tw3_im[2];
    y[2] = x0r + (w1r * x1r - w1i * x1i) + (w2r * x2r - w2i * x2i);
    y[3] = x0i + (w1r * x1i + w1i * x1r) + (w2r * x2i + w2i * x2r);
  }

  /* Y(2) = x0 + W3^2*x1 + W3*x2 */
  {
    const float32_t w1r = tw3_re[2], w1i = tw3_im[2];
    const float32_t w2r = tw3_re[1], w2i = tw3_im[1];
    y[4] = x0r + (w1r * x1r - w1i * x1i) + (w2r * x2r - w2i * x2i);
    y[5] = x0i + (w1r * x1i + w1i * x1r) + (w2r * x2i + w2i * x2r);
  }
}

/* Complex multiply: (a+ji)*(b+jc) -> out */
static inline void cmplx_mul(float32_t a, float32_t ai, float32_t b, float32_t bi,
                             float32_t *out_r, float32_t *out_i)
{
  *out_r = a * b - ai * bi;
  *out_i = a * bi + ai * b;
}

void fft_192_cplx_f32(const float32_t *x_384, float32_t *y_384,
                      const arm_cfft_instance_f32 *cfft_64)
{
  init_twiddles();

  /* Reshape: linear 192 complex -> 3x64 row-major (row i = samples i*64 .. i*64+63) */
  for (uint32_t k1 = 0; k1 < N1; k1++)
  {
    for (uint32_t n2 = 0; n2 < N2; n2++)
    {
      uint32_t n = k1 * N2 + n2;
      work_3x64[k1 * (N2 * 2) + n2 * 2 + 0] = x_384[n * 2 + 0];
      work_3x64[k1 * (N2 * 2) + n2 * 2 + 1] = x_384[n * 2 + 1];
    }
  }

  /* 64 x 3-point DFT on columns, then twiddle multiply (W3^(k1*n2)) */
  {
    float32_t col_in[6], col_out[6];
    for (uint32_t n2 = 0; n2 < N2; n2++)
    {
      /* Read column n2: 3 complex from the 3 rows */
      for (uint32_t k1 = 0; k1 < N1; k1++)
      {
        uint32_t off = k1 * (N2 * 2) + n2 * 2;
        col_in[k1 * 2 + 0] = work_3x64[off + 0];
        col_in[k1 * 2 + 1] = work_3x64[off + 1];
      }
      dft3_cplx_f32(col_in, col_out);
      /* Twiddle: multiply col_out[k1] by tw3[(k1*n2) % 3], write back */
      for (uint32_t k1 = 0; k1 < N1; k1++)
      {
        uint32_t tix = (k1 * n2) % 3;
        float32_t tr = tw3_re[tix], ti = tw3_im[tix];
        float32_t ar = col_out[k1 * 2 + 0], ai = col_out[k1 * 2 + 1];
        float32_t or_, oi;
        cmplx_mul(ar, ai, tr, ti, &or_, &oi);
        work_3x64[k1 * (N2 * 2) + n2 * 2 + 0] = or_;
        work_3x64[k1 * (N2 * 2) + n2 * 2 + 1] = oi;
      }
    }
  }

  /* 3 x 64-point complex FFT (each row) */
  for (uint32_t k1 = 0; k1 < N1; k1++)
    arm_cfft_f32(cfft_64, &work_3x64[k1 * (N2 * 2)], 0, 1);

  /* Unscramble: output bin k = k1 + k2*3 from row k1 bin k2 */
  for (uint32_t k = 0; k < N_FFT; k++)
  {
    uint32_t k1 = k % N1;
    uint32_t k2 = k / N1;
    y_384[k * 2 + 0] = work_3x64[k1 * (N2 * 2) + k2 * 2 + 0];
    y_384[k * 2 + 1] = work_3x64[k1 * (N2 * 2) + k2 * 2 + 1];
  }
}

void fft_192_real_f32(const float32_t *x_192, float32_t *y_384,
                      const arm_cfft_instance_f32 *cfft_64)
{
  /* Pack 192 real into 192 complex (imag = 0) then call complex FFT */
  __attribute__((aligned(32))) static float32_t x_cplx[FFT_192_OUTPUT_LEN];
  for (uint32_t i = 0; i < FFT_192_INPUT_LEN; i++)
  {
    x_cplx[i * 2 + 0] = x_192[i];
    x_cplx[i * 2 + 1] = 0.0f;
  }
  fft_192_cplx_f32(x_cplx, y_384, cfft_64);
}
