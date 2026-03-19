/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2026 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include "arm_math.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

/*
#define FFT_SIZE 256
#define ACTUAL_SAMPLES 185  // Number of real samples per chirp
*/
#define FFT_SIZE 192
#define ACTUAL_SAMPLES 192  // Number of real samples per chirp

#define NUM_CHIRPS 32       // Number of chirps to process
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
ADC_HandleTypeDef hadc1;
DMA_HandleTypeDef hdma_adc1;

I2C_HandleTypeDef hi2c4;

I2S_HandleTypeDef hi2s6;

LTDC_HandleTypeDef hltdc;

OSPI_HandleTypeDef hospi1;

RTC_HandleTypeDef hrtc;

SPI_HandleTypeDef hspi2;

UART_HandleTypeDef huart1;

SDRAM_HandleTypeDef hsdram1;

/* USER CODE BEGIN PV */
arm_cfft_instance_f32 cfft_instance;
arm_rfft_fast_instance_f32 rfft_64_instance;  // 64-point RFFT for 192-point (3×64) Range FFT

// OPTION B: FFT buffers in AXI SRAM for direct DMA-to-FFT flow (no copy overhead)
// With D-Cache enabled, performance is still excellent (~10-12K cycles vs 8.9K in DTCM)
// More realistic for continuous ADC streaming operation
__attribute__((aligned(32))) float32_t chirp_data[NUM_CHIRPS][ACTUAL_SAMPLES];
__attribute__((aligned(32))) float32_t fft_input[FFT_SIZE];
__attribute__((aligned(32))) float32_t fft_output[FFT_SIZE];
__attribute__((aligned(32))) float32_t cfft_buffer[NUM_CHIRPS * 2]; // Complex FFT working buffer (256 bytes)

// Large intermediate storage in regular RAM (AXI SRAM - fast but not as fast as DTCM)
__attribute__((aligned(32))) float32_t fft_outputs[NUM_CHIRPS][FFT_SIZE];       // 32 FFT results (32 KB)
__attribute__((aligned(32))) float32_t range_bins[FFT_SIZE/2][NUM_CHIRPS * 2];
__attribute__((aligned(32))) float32_t doppler_fft[FFT_SIZE/2][NUM_CHIRPS * 2];
__attribute__((aligned(32))) float32_t fft_magnitude[FFT_SIZE/2];                // 512 bytes

// Temporary arrays for optimized peak detection
__attribute__((aligned(32))) float32_t rd_map_magnitude_sq[FFT_SIZE/2][NUM_CHIRPS];  // Pre-computed |FFT|^2 for all cells (8 KB)
__attribute__((aligned(32))) float32_t rd_map_windowed_energy[FFT_SIZE/2][NUM_CHIRPS];  // Pre-computed 5×5 windowed energy (8 KB)
__attribute__((aligned(32))) float32_t rd_map_windowed_temp[FFT_SIZE/2][NUM_CHIRPS];  // Pre-computed 5×5 windowed energy (8 KB)

// ADC DMA buffer for continuous data acquisition (2 MSPS at 16-bit)
// Using circular buffer - DMA will continuously overwrite this
#define ADC_BUFFER_SIZE 8192  // 8K samples = ~4ms at 2 MSPS
__attribute__((aligned(32))) uint16_t adc_buffer[ADC_BUFFER_SIZE];  // 16 KB circular DMA buffer

// Hamming window coefficients for Range FFT (real FFT)
// Reduces spectral leakage by tapering signal at edges
__attribute__((aligned(32))) float32_t hamming_window_range[FFT_SIZE];

// Hamming window coefficients for Doppler FFT (complex FFT)
// Applied to 32 complex samples (across chirps for each range bin)
__attribute__((aligned(32))) float32_t hamming_window_doppler[NUM_CHIRPS];

// 192-point Range FFT (3×64) working buffers and twiddles
// - 3 complex buffers, each holds 64 complex samples (Re,Im interleaved) => 128 floats each
// - In DTCM for zero-wait access during the three arm_cfft_f32() calls (reduces Range 3xCFFT Total)
__attribute__((aligned(32), section(".dtcm_ram"))) float32_t range_fft_64_buf0[64 * 2];
__attribute__((aligned(32), section(".dtcm_ram"))) float32_t range_fft_64_buf1[64 * 2];
__attribute__((aligned(32), section(".dtcm_ram"))) float32_t range_fft_64_buf2[64 * 2];
// 3× 64 real buffers for RFFT input (demux); one packed temp for RFFT output before expand
__attribute__((aligned(32), section(".dtcm_ram"))) float32_t range_rfft_real0[64];
__attribute__((aligned(32), section(".dtcm_ram"))) float32_t range_rfft_real1[64];
__attribute__((aligned(32), section(".dtcm_ram"))) float32_t range_rfft_real2[64];
__attribute__((aligned(32), section(".dtcm_ram"))) float32_t range_rfft_packed[64];
// Twiddles for W192^k and W192^(2k) for k=0..(N/2), stored as complex (Re,Im)
// For N=192, this covers k=0..96 inclusive (97 values).
__attribute__((aligned(32))) float32_t range_twiddle_w1[(FFT_SIZE/2 + 1) * 2];
__attribute__((aligned(32))) float32_t range_twiddle_w2[(FFT_SIZE/2 + 1) * 2];

// Range FFT timing (overall per chirp and totals)
uint32_t range_fft_cycles;
float32_t range_fft_time_us;
uint32_t total_range_fft_cycles;    // Total cycles for all 32 range FFTs
float32_t total_range_fft_time_us;

// Detailed Range FFT timing (demux, 3×64 CFFTs, combine/pack)
uint32_t total_range_demux_cycles;      // Sum over all chirps
uint32_t total_range_cfft_cycles;       // Sum of three 64-pt CFFTs over all chirps
uint32_t total_range_combine_cycles;    // Sum of combine + pack over all chirps

uint32_t range_demux_cycles;            // Avg per chirp
uint32_t range_cfft_cycles;             // Avg per chirp
uint32_t range_combine_cycles;          // Avg per chirp

float32_t range_demux_time_us;          // Avg per chirp
float32_t range_cfft_time_us;           // Avg per chirp
float32_t range_combine_time_us;        // Avg per chirp

float32_t total_range_demux_time_us;    // Totals over all chirps
float32_t total_range_cfft_time_us;
float32_t total_range_combine_time_us;
uint32_t total_doppler_fft_cycles;  // Total cycles for all 64 Doppler FFTs (128-point FFT has 64 bins)
float32_t total_doppler_fft_time_us;      // Average Doppler FFT time in microseconds
uint32_t doppler_fft_cycles;        // Average cycles per Doppler FFT
float32_t doppler_fft_time_us;      // Average Doppler FFT time in microseconds

// Range windowing performance metrics
uint32_t total_range_window_cycles;       // Total cycles for all 32 range window applications
uint32_t range_window_cycles;             // Average cycles per range window application
float32_t range_window_time_us;           // Average range window time in microseconds
float32_t total_range_window_time_us;           // Average range window time in microseconds

// Doppler windowing performance metrics
uint32_t total_doppler_window_cycles;     // Total cycles for all 64 Doppler window applications
uint32_t doppler_window_cycles;           // Average cycles per Doppler window application
float32_t doppler_window_time_us;         // Average Doppler window time in microseconds
float32_t total_doppler_window_time_us;         // Average Doppler window time in microseconds

// Range-Doppler map peak detection results (5x5 windowed energy)
float32_t rd_map_peak_energy;             // Maximum windowed energy in Range-Doppler map
uint32_t rd_map_peak_range_bin;           // Range bin index of peak center (0-63)
uint32_t rd_map_peak_velocity_bin;        // Velocity bin index of peak center (0-31)
float32_t rd_map_peak_range_m;            // Actual range in meters
float32_t rd_map_peak_velocity_mps;       // Actual velocity in m/s

// Peak detection performance metrics
uint32_t peak_detection_cycles;           // Total cycles for peak detection
float32_t peak_detection_time_us;         // Peak detection time in microseconds

// Total time
uint32_t total_alg_cycles;           // Total cycles for peak detection
float32_t total_alg_time_us;         // Peak detection time in microseconds

// ADC DMA verification variables
uint16_t adc_value_before;          // ADC buffer value before DMA starts
uint16_t adc_value_after;           // ADC buffer value after DMA starts
uint32_t adc_changes;               // Count of buffer locations changed by DMA
/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
void PeriphCommonClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_DMA_Init(void);
static void MX_FMC_Init(void);
static void MX_I2C4_Init(void);
static void MX_I2S6_Init(void);
static void MX_LTDC_Init(void);
static void MX_OCTOSPI1_Init(void);
static void MX_RTC_Init(void);
static void MX_SPI2_Init(void);
static void MX_USART1_UART_Init(void);
static void MX_ADC1_Init(void);
/* USER CODE BEGIN PFP */
// static void MX_DMA_Init(void);   // Disabled - for ADC
// static void MX_ADC1_Init(void);  // Disabled - requires ADC HAL driver
/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

// ITM redirect for printf - sends output to SWV/ITM console
int _write(int file, char *ptr, int len)
{
  for (int i = 0; i < len; i++)
  {
    ITM_SendChar((*ptr++));
  }
  return len;
}

static inline void cmul_f32(float32_t a_re, float32_t a_im,
                            float32_t b_re, float32_t b_im,
                            float32_t *out_re, float32_t *out_im)
{
  *out_re = a_re * b_re - a_im * b_im;
  *out_im = a_re * b_im + a_im * b_re;
}

/* Twiddle + 3-point DFT for one k2; output bins written by caller (branch-free loops). */
static inline void range_fft_combine_one_k2(uint32_t k2,
                                            float32_t *out_y0_re, float32_t *out_y0_im,
                                            float32_t *out_y1_re, float32_t *out_y1_im)
{
  const float32_t w_re = -0.5f;
  const float32_t w_im = -0.8660254037844386f;
  const float32_t w2_re = -0.5f;
  const float32_t w2_im = 0.8660254037844386f;

  const float32_t x0_re = range_fft_64_buf0[2u * k2 + 0u];
  const float32_t x0_im = range_fft_64_buf0[2u * k2 + 1u];
  const float32_t x1_re0 = range_fft_64_buf1[2u * k2 + 0u];
  const float32_t x1_im0 = range_fft_64_buf1[2u * k2 + 1u];
  const float32_t x2_re0 = range_fft_64_buf2[2u * k2 + 0u];
  const float32_t x2_im0 = range_fft_64_buf2[2u * k2 + 1u];

  const float32_t tw1_re = range_twiddle_w1[2u * k2 + 0u];
  const float32_t tw1_im = range_twiddle_w1[2u * k2 + 1u];
  const float32_t tw2_re = range_twiddle_w2[2u * k2 + 0u];
  const float32_t tw2_im = range_twiddle_w2[2u * k2 + 1u];

  float32_t a1_re, a1_im, a2_re, a2_im;
  cmul_f32(x1_re0, x1_im0, tw1_re, tw1_im, &a1_re, &a1_im);
  cmul_f32(x2_re0, x2_im0, tw2_re, tw2_im, &a2_re, &a2_im);

  *out_y0_re = x0_re + a1_re + a2_re;
  *out_y0_im = x0_im + a1_im + a2_im;

  float32_t a1w_re, a1w_im, a1w2_re, a1w2_im;
  float32_t a2w_re, a2w_im, a2w2_re, a2w2_im;
  cmul_f32(a1_re, a1_im, w_re, w_im, &a1w_re, &a1w_im);
  cmul_f32(a1_re, a1_im, w2_re, w2_im, &a1w2_re, &a1w2_im);
  cmul_f32(a2_re, a2_im, w_re, w_im, &a2w_re, &a2w_im);
  cmul_f32(a2_re, a2_im, w2_re, w2_im, &a2w2_re, &a2w2_im);

  *out_y1_re = x0_re + a1w_re + a2w2_re;
  *out_y1_im = x0_im + a1w_im + a2w2_im;
}

// Expand 64-point RFFT packed output to 128-float complex (Re,Im interleaved) for combine.
// packed: CMSIS layout [0]=DC, [1]=Nyquist, [2*k],[2*k+1]=Re(X[k]),Im(X[k]) for k=1..31.
// complex_buf: [2*k2]=Re(S[k2]), [2*k2+1]=Im(S[k2]) for k2=0..63; bins 33..63 from conjugate symmetry.
static void range_rfft_expand_64_to_128(const float32_t *packed, float32_t *complex_buf)
{
  complex_buf[0u] = packed[0u];
  complex_buf[1u] = 0.0f;
  for (uint32_t k2 = 1u; k2 < 32u; k2++)
  {
    complex_buf[2u * k2 + 0u] = packed[2u * k2 + 0u];
    complex_buf[2u * k2 + 1u] = packed[2u * k2 + 1u];
  }
  complex_buf[64u] = packed[1u];
  complex_buf[65u] = 0.0f;
  for (uint32_t k2 = 33u; k2 < 64u; k2++)
  {
    const uint32_t k2_conj = 64u - k2;
    complex_buf[2u * k2 + 0u] = complex_buf[2u * k2_conj + 0u];
    complex_buf[2u * k2 + 1u] = -complex_buf[2u * k2_conj + 1u];
  }
}

// 192-point real-input FFT: Cooley–Tukey N=3×64 with time index n = n1 + 3*n2
// (polyphase-by-3 demux). Three 64-pt RFFTs then expand + twiddle + 3-point DFT.
// X[k2 + 64*r] = sum_{n1=0}^{2} F[n1,k2] * W_192^{(k2 + 64*r)*n1}, r=0,1,2.
// Output packed like CMSIS `arm_rfft_fast_f32`: DC, Nyquist, then Re/Im pairs.
static void range_fft_192_real_packed(const float32_t *time_in, float32_t *packed_out)
{
  uint32_t t0 = DWT->CYCCNT;

  // Demux: seq_r[n2] = x[n1 + 3*n2] with n1=r, n2=0..63 into three real length-64 buffers
  for (uint32_t n2 = 0u; n2 < 64u; n2++)
  {
    range_rfft_real0[n2] = time_in[3u * n2 + 0u];
    range_rfft_real1[n2] = time_in[3u * n2 + 1u];
    range_rfft_real2[n2] = time_in[3u * n2 + 2u];
  }

  uint32_t t1 = DWT->CYCCNT;

  arm_rfft_fast_f32(&rfft_64_instance, range_rfft_real0, range_rfft_packed, 0);
  range_rfft_expand_64_to_128(range_rfft_packed, range_fft_64_buf0);
  arm_rfft_fast_f32(&rfft_64_instance, range_rfft_real1, range_rfft_packed, 0);
  range_rfft_expand_64_to_128(range_rfft_packed, range_fft_64_buf1);
  arm_rfft_fast_f32(&rfft_64_instance, range_rfft_real2, range_rfft_packed, 0);
  range_rfft_expand_64_to_128(range_rfft_packed, range_fft_64_buf2);

  uint32_t t2 = DWT->CYCCNT;

  float32_t y0_re, y0_im, y1_re, y1_im;

  range_fft_combine_one_k2(0u, &y0_re, &y0_im, &y1_re, &y1_im);
  packed_out[0] = y0_re;
  packed_out[128u] = y1_re;
  packed_out[129u] = y1_im;

  for (uint32_t k2 = 1u; k2 < 32u; k2++)
  {
    range_fft_combine_one_k2(k2, &y0_re, &y0_im, &y1_re, &y1_im);
    packed_out[2u * k2 + 0u] = y0_re;
    packed_out[2u * k2 + 1u] = y0_im;
    const uint32_t kb = k2 + 64u;
    packed_out[2u * kb + 0u] = y1_re;
    packed_out[2u * kb + 1u] = y1_im;
  }

  range_fft_combine_one_k2(32u, &y0_re, &y0_im, &y1_re, &y1_im);
  packed_out[64u] = y0_re;
  packed_out[65u] = y0_im;
  packed_out[1] = y1_re; /* Re{X[96]} Nyquist */

  for (uint32_t k2 = 33u; k2 < 64u; k2++)
  {
    range_fft_combine_one_k2(k2, &y0_re, &y0_im, &y1_re, &y1_im);
    packed_out[2u * k2 + 0u] = y0_re;
    packed_out[2u * k2 + 1u] = y0_im;
  }

  uint32_t t3 = DWT->CYCCNT;

  total_range_demux_cycles   += (t1 - t0);
  total_range_cfft_cycles    += (t2 - t1);
  total_range_combine_cycles += (t3 - t2);
}

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{

  /* USER CODE BEGIN 1 */
	volatile static uint32_t icache_enabled = 0;
	volatile static uint32_t dcache_enabled = 0;
  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */
  // Enable CPU I-Cache and D-Cache for maximum performance
  SCB_EnableICache();
  SCB_EnableDCache();

  // Verify caches are enabled (for debugging)
  icache_enabled = (SCB->CCR & SCB_CCR_IC_Msk) ? 1 : 0;
  dcache_enabled = (SCB->CCR & SCB_CCR_DC_Msk) ? 1 : 0;

  // Verify arrays are in DTCM (should be 0x20000000-0x2001FFFF range)
  uint32_t input_addr = (uint32_t)&fft_input[0];
  uint32_t output_addr = (uint32_t)&fft_output[0];
  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* Configure the peripherals common clocks */
  PeriphCommonClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_DMA_Init();
  MX_FMC_Init();
  MX_I2C4_Init();
  MX_I2S6_Init();
  MX_LTDC_Init();
  MX_OCTOSPI1_Init();
  MX_RTC_Init();
  MX_SPI2_Init();
  MX_USART1_UART_Init();
  MX_ADC1_Init();
  /* USER CODE BEGIN 2 */

  printf("icache_enabled = %d. dcache_enabled=%d\r\n", (int)icache_enabled, (int) dcache_enabled);
  printf("input_addr = %x. output_addr=%x\r\n", (unsigned int)input_addr, (unsigned int) output_addr);


  // Enable DWT cycle counter for time measurement
  CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
  DWT->CYCCNT = 0;
  DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;

  // Initialize FFT instances
  arm_cfft_init_f32(&cfft_instance, NUM_CHIRPS);          // 32-point complex FFT for Doppler
  arm_rfft_fast_init_f32(&rfft_64_instance, 64);          // 64-point real FFT for 192-point (3×64) Range FFT

  // Generate Hamming window coefficients for Range FFT (FFT_SIZE points)
  // Hamming window: w(n) = 0.54 - 0.46 * cos(2*pi*n/(N-1))
  // Reduces spectral leakage by tapering signal at edges
  for (uint32_t n = 0; n < FFT_SIZE; n++)
  {
    hamming_window_range[n] = 0.54f - 0.46f * arm_cos_f32(2.0f * PI * (float32_t)n / (float32_t)(FFT_SIZE - 1));
  }

  // Generate Hamming window coefficients for Doppler FFT (32 points)
  // Applied across chirps for velocity resolution
  for (uint32_t n = 0; n < NUM_CHIRPS; n++)
  {
    hamming_window_doppler[n] = 0.54f - 0.46f * arm_cos_f32(2.0f * PI * (float32_t)n / (float32_t)(NUM_CHIRPS - 1));
  }

  // Precompute twiddles for the 192-point (3×64) Range FFT combine step.
  // W192^k = exp(-j*2*pi*k/192) and W192^(2k) = exp(-j*4*pi*k/192)
  for (uint32_t k = 0; k <= (FFT_SIZE/2); k++)
  {
    float32_t a1 = -2.0f * PI * (float32_t)k / (float32_t)FFT_SIZE;
    float32_t a2 = -4.0f * PI * (float32_t)k / (float32_t)FFT_SIZE;

    range_twiddle_w1[2u * k + 0u] = arm_cos_f32(a1);
    range_twiddle_w1[2u * k + 1u] = arm_sin_f32(a1);

    range_twiddle_w2[2u * k + 0u] = arm_cos_f32(a2);
    range_twiddle_w2[2u * k + 1u] = arm_sin_f32(a2);
  }

  // ===================================================================
  // Start ADC continuous acquisition with DMA (2 MSPS at 16-bit)
  // ===================================================================
  // This continuously samples ADC and DMAs to circular buffer
  // Simulates real-world radar ADC acquisition during FFT processing
  // DMA will compete for memory bus bandwidth with FFT
  // This tests realistic performance under memory bus contention

  // Clear buffer to verify DMA is writing
  for (uint32_t i = 0; i < ADC_BUFFER_SIZE; i++)
  {
    adc_buffer[i] = 0xDEAD;  // Marker pattern
  }

  adc_value_before = adc_buffer[100];  // Sample before starting DMA
  adc_changes = 0;  // Count how many buffer locations changed

  HAL_ADC_Start_DMA(&hadc1, (uint32_t*)adc_buffer, ADC_BUFFER_SIZE);

  // Small delay to allow DMA to start transferring
  for (volatile uint32_t i = 0; i < 100000; i++);

  adc_value_after = adc_buffer[100];   // Sample after DMA starts

  // Count how many values changed from the marker pattern
  for (uint32_t i = 0; i < ADC_BUFFER_SIZE; i++)
  {
    if (adc_buffer[i] != 0xDEAD)
      adc_changes++;
  }

  // ===================================================================
  // Acquire ADC data for all chirps
  // ===================================================================
  // In real radar: DMA continuously fills buffer, we extract chirps
  // For testing: Convert ADC buffer samples to float32_t chirp data
  // ADC is 16-bit: 0-65535, normalize to ±1.0 range (centered at 32768)

  float32_t sampling_freq = 2000000.0f;  // 2 MSPS actual ADC rate

  // Extract chirps from circular ADC buffer
  // Each chirp is ACTUAL_SAMPLES (185) consecutive samples
  uint32_t adc_index = 0;  // Starting position in ADC buffer

#define USE_TEST_SIGNAL 1  // Set to 0 to use real ADC input

#if USE_TEST_SIGNAL
  // ===================================================================
  // REALISTIC FMCW RADAR SIGNAL SIMULATION
  // ===================================================================
  // Simulates a radar target with specific range and velocity
  // This creates proper Range-Doppler response for system validation

  // ===== Radar Parameters =====
  const float32_t c = 3.0e8f;                          // Speed of light (m/s)
  const float32_t radar_freq_hz = 77.0e9f;             // Radar center frequency: 77 GHz (automotive radar)
  const float32_t chirp_bandwidth_hz = 150.0e6f;       // Chirp bandwidth: 150 MHz
  const float32_t chirp_duration_s = (float32_t)ACTUAL_SAMPLES / sampling_freq;  // 64 µs for 128 samples at 2 MSPS
  const float32_t chirp_repetition_interval_s = 100.0e-6f;  // 100 µs between chirps (PRI)

  // ===== Simulated Target Parameters =====
  const float32_t target_range_m = 50.0f;              // Target at 50 meters
  const float32_t target_velocity_mps = 8.0f;          // Target velocity: 8 m/s (28.8 km/h) approaching (within unambiguous range)

  // Calculate beat frequency from target range
  // For FMCW: f_beat = 2 * R * B / (c * T_chirp)
  float32_t beat_frequency = 2.0f * target_range_m * chirp_bandwidth_hz / (c * chirp_duration_s);

  // Calculate Doppler frequency from target velocity
  // f_doppler = 2 * v * f_c / c
  float32_t doppler_frequency = 2.0f * target_velocity_mps * radar_freq_hz / c;

  // Generate FMCW radar signal for all chirps
  for (uint32_t chirp = 0; chirp < NUM_CHIRPS; chirp++)
  {
    // Slow-time (between chirps) for Doppler phase progression
    float32_t t_slow = (float32_t)chirp * chirp_repetition_interval_s;

    // Doppler phase changes between chirps
    float32_t doppler_phase = 2.0f * PI * doppler_frequency * t_slow;

    for (uint32_t i = 0; i < ACTUAL_SAMPLES; i++)
    {
      // Fast-time (within chirp) for range beat frequency
      float32_t t_fast = (float32_t)i / sampling_freq;

      // Combined radar signal: beat frequency (range) + Doppler phase (velocity)
      // Using cosine for I-channel (in real radar, you'd also have Q-channel)
      chirp_data[chirp][i] = arm_cos_f32(2.0f * PI * beat_frequency * t_fast + doppler_phase);

      adc_index++;
      if (adc_index >= ADC_BUFFER_SIZE)
        adc_index = 0;
    }
  }

  // Calculate expected bin locations for verification
  float32_t range_resolution = c / (2.0f * chirp_bandwidth_hz);  // Meters per bin
  float32_t expected_range_bin = target_range_m / range_resolution;

  float32_t velocity_resolution = c / (2.0f * radar_freq_hz * NUM_CHIRPS * chirp_repetition_interval_s);  // m/s per bin
  float32_t expected_velocity_bin = target_velocity_mps / velocity_resolution;  // Bin in complex FFT (0-15 positive, 16-31 negative)

  // Print expected results (for verification)
  printf("\r\n--- Simulated Radar Target ---\r\n");
  printf("Target Range:     %.1f m\r\n", target_range_m);
  printf("Target Velocity:  %.1f m/s (%.1f km/h)\r\n", target_velocity_mps, target_velocity_mps * 3.6f);
  printf("Beat Frequency:   %.2f Hz\r\n", beat_frequency);
  printf("Doppler Freq:     %.2f Hz\r\n", doppler_frequency);
  printf("Expected Range Bin:    %.1f\r\n", expected_range_bin);
  printf("Expected Velocity Bin: %.1f\r\n", expected_velocity_bin);
  printf("Range Resolution: %.2f m/bin\r\n", range_resolution);
  printf("Velocity Resolution: %.2f m/s/bin\r\n", velocity_resolution);
  printf("==============================\r\n\r\n");
#else
  // ===================================================================
  // OPTION 2: Use real ADC input (requires signal connected to ADC pin)
  // ===================================================================
  for (uint32_t chirp = 0; chirp < NUM_CHIRPS; chirp++)
  {
    for (uint32_t i = 0; i < ACTUAL_SAMPLES; i++)
    {
      // Convert ADC sample (uint16_t 0-65535) to float32_t normalized to ±1.0
      // Center at 32768 (mid-scale) and scale to ±1.0
      int32_t adc_signed = (int32_t)adc_buffer[adc_index] - 32768;
      chirp_data[chirp][i] = (float32_t)adc_signed / 32768.0f;

      // Advance to next sample, wrap around circular buffer
      adc_index++;
      if (adc_index >= ADC_BUFFER_SIZE)
        adc_index = 0;
    }
  }
#endif

  // === STEP 1: Range FFT (32 real FFTs on chirps) ===

  // Warm-up run to ensure code and data are in cache
  for (uint32_t i = 0; i < ACTUAL_SAMPLES; i++)
  {
    fft_input[i] = chirp_data[0][i];
  }
  for (uint32_t i = ACTUAL_SAMPLES; i < FFT_SIZE; i++)
  {
    fft_input[i] = 0.0f;
  }
  // Apply Hamming window to warm-up data
  for (uint32_t i = 0; i < FFT_SIZE; i++)
  {
    fft_input[i] *= hamming_window_range[i];
  }
  range_fft_192_real_packed(fft_input, fft_output);

  // Process all 32 chirps and measure windowing and FFT computation time separately
  total_range_fft_cycles = 0;
  total_range_window_cycles = 0;
  total_range_demux_cycles = 0;
  total_range_cfft_cycles = 0;
  total_range_combine_cycles = 0;

  for (uint32_t chirp = 0; chirp < NUM_CHIRPS; chirp++)
  {
    // Copy chirp data to FFT input buffer - NOT TIMED
    for (uint32_t i = 0; i < ACTUAL_SAMPLES; i++)
    {
      fft_input[i] = chirp_data[chirp][i];
    }

    // Zero-pad the remaining samples (if ACTUAL_SAMPLES < FFT_SIZE) - NOT TIMED
    for (uint32_t i = ACTUAL_SAMPLES; i < FFT_SIZE; i++)
    {
      fft_input[i] = 0.0f;
    }

    // Apply Range Hamming window and measure time
    uint32_t window_start_cycles = DWT->CYCCNT;
    for (uint32_t i = 0; i < FFT_SIZE; i++)
    {
      fft_input[i] *= hamming_window_range[i];
    }
    uint32_t window_end_cycles = DWT->CYCCNT;
    total_range_window_cycles += (window_end_cycles - window_start_cycles);

    // Measure ONLY the FFT computation
    uint32_t start_cycles = DWT->CYCCNT;
    range_fft_192_real_packed(fft_input, fft_output);
    uint32_t end_cycles = DWT->CYCCNT;
    total_range_fft_cycles += (end_cycles - start_cycles);

    // Store FFT result - NOT TIMED
    for (uint32_t i = 0; i < FFT_SIZE; i++)
    {
      fft_outputs[chirp][i] = fft_output[i];
    }
  }

  // Calculate average range FFT execution time per chirp
  range_fft_cycles = total_range_fft_cycles / NUM_CHIRPS;
  range_fft_time_us = (float32_t)range_fft_cycles / 280.0f;
  total_range_fft_time_us = (float32_t)total_range_fft_cycles / 280.0f;

  // Calculate average range windowing execution time per chirp
  range_window_cycles = total_range_window_cycles / NUM_CHIRPS;
  range_window_time_us = (float32_t)range_window_cycles / 280.0f;
  total_range_window_time_us = (float32_t)total_range_window_cycles / 280.0f;

  // Calculate detailed Range FFT breakdown per chirp and totals
  range_demux_cycles     = total_range_demux_cycles / NUM_CHIRPS;
  range_cfft_cycles      = total_range_cfft_cycles / NUM_CHIRPS;
  range_combine_cycles   = total_range_combine_cycles / NUM_CHIRPS;

  range_demux_time_us    = (float32_t)range_demux_cycles / 280.0f;
  range_cfft_time_us     = (float32_t)range_cfft_cycles / 280.0f;
  range_combine_time_us  = (float32_t)range_combine_cycles / 280.0f;

  total_range_demux_time_us   = (float32_t)total_range_demux_cycles / 280.0f;
  total_range_cfft_time_us    = (float32_t)total_range_cfft_cycles / 280.0f;
  total_range_combine_time_us = (float32_t)total_range_combine_cycles / 280.0f;

  // Reorganize FFT outputs from [chirp][bin] to [bin][chirp] format
  // This creates a range-Doppler matrix suitable for 2D FFT processing
  for (uint32_t bin = 0; bin < FFT_SIZE/2; bin++)
  {
    for (uint32_t chirp = 0; chirp < NUM_CHIRPS; chirp++)
    {
      if (bin == 0)
      {
        // DC component (bin 0): real only, imaginary = 0
        range_bins[0][chirp * 2 + 0] = fft_outputs[chirp][0];  // Real
        range_bins[0][chirp * 2 + 1] = 0.0f;                    // Imag = 0
      }
      else
      {
        // Bins 1-127: complex pairs starting at position [2]
        range_bins[bin][chirp * 2 + 0] = fft_outputs[chirp][bin * 2];      // Real
        range_bins[bin][chirp * 2 + 1] = fft_outputs[chirp][bin * 2 + 1];  // Imag
      }
    }
  }

  // range_bins[128][64] now contains reorganized data:
  // range_bins[bin][0,1] = complex sample from chirp 0, bin N (real, imag)
  // range_bins[bin][2,3] = complex sample from chirp 1, bin N (real, imag)
  // ...
  // range_bins[bin][62,63] = complex sample from chirp 31, bin N (real, imag)

  // === STEP 3: Doppler FFT (128 complex FFTs across chirps) ===

  // Warm-up run for complex FFT with windowing
  for (uint32_t i = 0; i < NUM_CHIRPS * 2; i++)
  {
    cfft_buffer[i] = range_bins[0][i];
  }
  // Apply Doppler Hamming window to warm-up data
  // For complex data: multiply both real and imaginary parts by same coefficient
  for (uint32_t i = 0; i < NUM_CHIRPS; i++)
  {
    cfft_buffer[i * 2 + 0] *= hamming_window_doppler[i];  // Real part
    cfft_buffer[i * 2 + 1] *= hamming_window_doppler[i];  // Imaginary part
  }
  arm_cfft_f32(&cfft_instance, cfft_buffer, 0, 1);

  // Process all 64 range bins and measure windowing and FFT computation time separately
  total_doppler_fft_cycles = 0;
  total_doppler_window_cycles = 0;

  for (uint32_t bin = 0; bin < FFT_SIZE/2; bin++)
  {
    // Copy range bin data to working buffer - NOT TIMED
    for (uint32_t i = 0; i < NUM_CHIRPS * 2; i++)
    {
      cfft_buffer[i] = range_bins[bin][i];
    }

    // Apply Doppler Hamming window and measure time
    // For complex data: multiply both real and imaginary parts by same coefficient
    uint32_t window_start_cycles = DWT->CYCCNT;
    for (uint32_t i = 0; i < NUM_CHIRPS; i++)
    {
      cfft_buffer[i * 2 + 0] *= hamming_window_doppler[i];  // Real part
      cfft_buffer[i * 2 + 1] *= hamming_window_doppler[i];  // Imaginary part
    }
    uint32_t window_end_cycles = DWT->CYCCNT;
    total_doppler_window_cycles += (window_end_cycles - window_start_cycles);

    // Measure ONLY the complex FFT computation
    uint32_t start_cycles = DWT->CYCCNT;
    arm_cfft_f32(&cfft_instance, cfft_buffer, 0, 1);  // 0 = forward FFT, 1 = bit reverse
    uint32_t end_cycles = DWT->CYCCNT;
    total_doppler_fft_cycles += (end_cycles - start_cycles);

    // Store Doppler FFT result - NOT TIMED
    for (uint32_t i = 0; i < NUM_CHIRPS * 2; i++)
    {
      doppler_fft[bin][i] = cfft_buffer[i];
    }
  }

  // Calculate average Doppler FFT execution time per bin
  doppler_fft_cycles = total_doppler_fft_cycles / (FFT_SIZE/2);
  doppler_fft_time_us = (float32_t)doppler_fft_cycles / 280.0f;
  total_doppler_fft_time_us = (float32_t)total_doppler_fft_cycles / 280.0f;

  // Calculate average Doppler windowing execution time per bin
  doppler_window_cycles = total_doppler_window_cycles / (FFT_SIZE/2);
  doppler_window_time_us = (float32_t)doppler_window_cycles / 280.0f;
  total_doppler_window_time_us = (float32_t)total_doppler_window_cycles / 280.0f;

  // Example: Calculate magnitude spectrum for first chirp
  arm_cmplx_mag_f32(fft_outputs[0], fft_magnitude, FFT_SIZE/2);

  // Find peak frequency in first chirp
  float32_t max_value;
  uint32_t max_index;
  arm_max_f32(fft_magnitude, FFT_SIZE/2, &max_value, &max_index);

  // Calculate frequency of peak (in Hz)
  float32_t peak_frequency = (float32_t)max_index * sampling_freq / (float32_t)FFT_SIZE;

  // =============================================================================
  // Find Peak in Range-Doppler Map (5×5 Windowed Energy - HIGHLY OPTIMIZED)
  // =============================================================================
  // Separable filtering approach for 5×5 window:
  // 1. Pre-compute magnitude² for all 2,048 cells (no sqrt!, unrolled by 4)
  // 2a. Horizontal 1×5 filtering using sliding window (2 ops per cell after init)
  // 2b. Vertical 5×1 filtering using sliding window (2 ops per cell after init)
  // 3. Simple max search over 1,680 windowed energy values
  //
  // Complexity: 8 operations per cell vs 24 in naive approach (3× reduction!)
  // Memory access: Sequential, cache-friendly patterns with restrict pointers

  uint32_t peak_start_cycles = DWT->CYCCNT;

  // STEP 1: Pre-compute magnitude squared for all Range-Doppler map cells
  // Optimized with restrict pointers for SIMD vectorization
  for (uint32_t r = 0; r < FFT_SIZE/2; r++)
  {
    // Use restrict to tell compiler arrays don't overlap (enables vectorization)
    float32_t * __restrict pSrc = doppler_fft[r];           // [Re0, Im0, Re1, Im1, ...]
    float32_t * __restrict pDst = rd_map_magnitude_sq[r];

    // Unroll by 4 for better SIMD utilization
    uint32_t v;
    for (v = 0; v < (NUM_CHIRPS & ~3); v += 4)
    {
      // Process 4 complex values at once
      float32_t r0 = pSrc[0], i0 = pSrc[1];
      float32_t r1 = pSrc[2], i1 = pSrc[3];
      float32_t r2 = pSrc[4], i2 = pSrc[5];
      float32_t r3 = pSrc[6], i3 = pSrc[7];

      pDst[0] = r0*r0 + i0*i0;
      pDst[1] = r1*r1 + i1*i1;
      pDst[2] = r2*r2 + i2*i2;
      pDst[3] = r3*r3 + i3*i3;

      pSrc += 8;
      pDst += 4;
    }

    // Handle remaining elements (NUM_CHIRPS = 32, so this won't execute)
    for (; v < NUM_CHIRPS; v++)
    {
      float32_t real = *pSrc++;
      float32_t imag = *pSrc++;
      *pDst++ = real*real + imag*imag;
    }
  }

  // STEP 2: Separable 5×5 windowed energy computation (8 adds per cell vs 24)
  // Step 2a: Horizontal 1×5 filtering
  // Step 2b: Vertical 5×1 filtering
  #define WINDOW_SIZE 5
  #define WINDOW_HALF (WINDOW_SIZE / 2)  // = 2

  #if WINDOW_SIZE != 5
  #error change sum formula loop code in step 2a and step 2b

  #endif
  // STEP 2a: Horizontal pass - compute 1×5 sum for each row
  // Process FULL rows to maximize inner loop iterations (better for CPU pipeline)
  uint32_t row_start_cycles = DWT->CYCCNT;
  for (uint32_t r = 0; r < FFT_SIZE/2; r++)
  {
    // Use restrict pointers to help compiler optimize
    float32_t * __restrict pSrc = rd_map_magnitude_sq[r];
    float32_t * __restrict pDst = rd_map_windowed_temp[r];

    // Process valid columns using sliding window
    // Initial window (columns 0-4)
    float32_t sum = pSrc[0] + pSrc[1] + pSrc[2] + pSrc[3] + pSrc[4];
    pDst[WINDOW_HALF] = sum;

    // Slide window across row (columns 3-31)
    for (uint32_t v = WINDOW_HALF + 1; v < (NUM_CHIRPS - WINDOW_HALF); v++)
    {
      sum = sum - pSrc[v - WINDOW_HALF - 1] + pSrc[v + WINDOW_HALF];
      pDst[v] = sum;
    }
  }
  uint32_t row_end_cycles = DWT->CYCCNT;
  uint32_t row_cycles = row_end_cycles - row_start_cycles;

  // STEP 2b: Vertical pass - compute 5×1 sum for each column
  // Process FULL columns
  uint32_t col_start_cycles = DWT->CYCCNT;
  for (uint32_t v = WINDOW_HALF; v < (NUM_CHIRPS - WINDOW_HALF); v++)
  {
    // Initial window (rows 0-4) for this column
    float32_t sum = rd_map_windowed_temp[0][v] +
                    rd_map_windowed_temp[1][v] +
                    rd_map_windowed_temp[2][v] +
                    rd_map_windowed_temp[3][v] +
                    rd_map_windowed_temp[4][v];
    rd_map_windowed_energy[WINDOW_HALF][v] = sum;

    // Slide window down column (rows 3-61)
    for (uint32_t r = WINDOW_HALF + 1; r < (FFT_SIZE/2 - WINDOW_HALF); r++)
    {
      sum = sum - rd_map_windowed_temp[r - WINDOW_HALF - 1][v] + rd_map_windowed_temp[r + WINDOW_HALF][v];
      rd_map_windowed_energy[r][v] = sum;
    }
  }

  uint32_t col_end_cycles = DWT->CYCCNT;
  uint32_t col_cycles = col_end_cycles - col_start_cycles;

  // STEP 3: Find maximum in windowed energy matrix
  // Optimized linear search with minimal branches
  rd_map_peak_energy = 0.0f;
  rd_map_peak_range_bin = WINDOW_HALF;
  rd_map_peak_velocity_bin = WINDOW_HALF;

  for (uint32_t r = WINDOW_HALF; r < (FFT_SIZE/2 - WINDOW_HALF); r++)
  {
    float32_t * __restrict pRow = rd_map_windowed_energy[r];

    for (uint32_t v = WINDOW_HALF; v < (NUM_CHIRPS - WINDOW_HALF); v++)
    {
      float32_t energy = pRow[v];
      // Branchless comparison using conditional move (compiler optimization)
      if (energy > rd_map_peak_energy)
      {
        rd_map_peak_energy = energy;
        rd_map_peak_range_bin = r;
        rd_map_peak_velocity_bin = v;
      }
    }
  }

  uint32_t peak_end_cycles = DWT->CYCCNT;
  peak_detection_cycles = peak_end_cycles - peak_start_cycles;
  peak_detection_time_us = (float32_t)peak_detection_cycles / 280.0f;

  total_alg_cycles = total_range_window_cycles + total_range_fft_cycles + total_doppler_window_cycles + total_doppler_fft_cycles + peak_detection_cycles;
  total_alg_time_us = (float32_t)total_alg_cycles / 280.0f;



  // Calculate actual physical values from bin indices using FMCW radar parameters
  // Must match the parameters used in signal generation (77 GHz, 150 MHz BW, 100 µs PRI)

  // Range calculation: Range resolution = c / (2 * Bandwidth)
  // With 150 MHz bandwidth: resolution = 3e8 / (2 * 150e6) = 1.0 meter/bin
  rd_map_peak_range_m = (float32_t)rd_map_peak_range_bin * (3.0e8f / (2.0f * 150.0e6f));

  // Velocity calculation: Velocity resolution = c / (2 * fc * N * PRI)
  // With 77 GHz, 32 chirps, 100 µs PRI: resolution ≈ 0.609 m/s per bin
  // Complex FFT: bin 0 = DC (0 m/s), bins 1-15 = positive freq, bins 16-31 = negative freq (wrap around)
  float32_t vel_res = 3.0e8f / (2.0f * 77.0e9f * (float32_t)NUM_CHIRPS * 100.0e-6f);
  if (rd_map_peak_velocity_bin < NUM_CHIRPS/2) {
    rd_map_peak_velocity_mps = (float32_t)rd_map_peak_velocity_bin * vel_res;  // Positive velocity
  } else {
    rd_map_peak_velocity_mps = ((float32_t)rd_map_peak_velocity_bin - (float32_t)NUM_CHIRPS) * vel_res;  // Negative velocity (wrapped)
  }

  // =============================================================================
  // Print Results to ITM Console (does NOT affect timing - all measurements complete)
  // =============================================================================

  printf("\r\n========================================\r\n");
  printf("RADAR FFT PROCESSING RESULTS\r\n");
  printf("========================================\r\n");

  printf("\r\n--- Timing Performance @ 280 MHz ---\r\n");
  printf("Range Window:   %6lu cycles (%7.2f us) [avg per chirp]\r\n", range_window_cycles, range_window_time_us);
  printf("Range FFT:      %6lu cycles (%7.2f us) [avg per chirp]\r\n", range_fft_cycles, range_fft_time_us);
  printf("  Range Demux:  %6lu cycles (%7.2f us) [avg per chirp]\r\n", range_demux_cycles, range_demux_time_us);
  printf("  Range 3xCFFT: %6lu cycles (%7.2f us) [avg per chirp]\r\n", range_cfft_cycles, range_cfft_time_us);
  printf("  Range Comb.:  %6lu cycles (%7.2f us) [avg per chirp]\r\n", range_combine_cycles, range_combine_time_us);
  printf("Doppler Window: %6lu cycles (%7.2f us) [avg per bin]\r\n", doppler_window_cycles, doppler_window_time_us);
  printf("Doppler FFT:    %6lu cycles (%7.2f us) [avg per bin]\r\n", doppler_fft_cycles, doppler_fft_time_us);

  printf("Range Window Total:   %6lu cycles (%7.2f us)\r\n", total_range_window_cycles, total_range_window_time_us);
  printf("Range FFT Total:      %6lu cycles (%7.2f us)\r\n", total_range_fft_cycles, total_range_fft_time_us);
  printf("  Range Demux Total:  %6lu cycles (%7.2f us)\r\n", total_range_demux_cycles, total_range_demux_time_us);
  printf("  Range 3xCFFT Total: %6lu cycles (%7.2f us)\r\n", total_range_cfft_cycles, total_range_cfft_time_us);
  printf("  Range Comb. Total:  %6lu cycles (%7.2f us)\r\n", total_range_combine_cycles, total_range_combine_time_us);
  printf("Doppler Window Total: %6lu cycles (%7.2f us)\r\n", total_doppler_window_cycles, total_doppler_window_time_us);
  printf("Doppler FFT Total:    %6lu cycles (%7.2f us)\r\n", total_doppler_fft_cycles, total_doppler_fft_time_us);
  printf("Peak Detection: %6lu cycles (%7.2f us)\r\n", peak_detection_cycles, peak_detection_time_us);
  printf("optimization: row %6lu cycles col %6lu cycles\r\n", row_cycles, col_cycles);
  printf("---\r\n");
  printf("TOTAL TIME:     %6lu cycles (%7.2f us = %.3f ms)\r\n",
         total_alg_cycles, total_alg_time_us, total_alg_time_us / 1000.0f);

  printf("\r\n--- Peak Detection Results ---\r\n");
  printf("Peak Energy:      %.2f\r\n", rd_map_peak_energy);
  printf("Range Bin:        %lu (%.2f m)\r\n", rd_map_peak_range_bin, rd_map_peak_range_m);
  printf("Velocity Bin:     %lu (%.2f m/s)\r\n", rd_map_peak_velocity_bin, rd_map_peak_velocity_mps);

  printf("\r\n--- Configuration ---\r\n");
  printf("FFT Size:         %d points\r\n", FFT_SIZE);
  printf("Samples/chirp:    %d\r\n", ACTUAL_SAMPLES);
  printf("Number of chirps: %d\r\n", NUM_CHIRPS);
  printf("ADC changed:      %lu samples\r\n", adc_changes);
  printf("Peak frequency:   %.2f Hz\r\n", peak_frequency);

  printf("========================================\r\n\r\n");

  // =============================================================================
  // PROCESSING COMPLETE - Set a breakpoint here to inspect results
  // =============================================================================
  //
  // SYSTEM ARCHITECTURE:
  // - ADC continuously samples at 2 MSPS (16-bit) via DMA to AXI SRAM
  // - DMA writes to adc_buffer[] in AXI SRAM (circular buffer)
  // - CPU extracts chirps and converts uint16_t → float32_t
  // - Hamming windows applied to both Range and Doppler FFTs
  //   * Range: 128-point window reduces spectral leakage in range dimension
  //   * Doppler: 32-point window reduces spectral leakage in velocity dimension
  // - FFT runs directly from AXI SRAM (no copy to DTCM)
  // - D-Cache enabled for good performance despite AXI SRAM access
  // - DMA and FFT share AXI bus → realistic memory contention
  //
  // EXPECTED PERFORMANCE (128-point Range FFT, 32-point Doppler FFT):
  // - Range Hamming window: ~100-300 cycles (128 float multiplications)
  // - Range FFT: ~4-6K cycles (128-point real FFT)
  // - Doppler Hamming window: ~50-100 cycles (64 float multiplications for complex data)
  // - Doppler FFT: ~2-3K cycles per 32-point complex FFT
  // - Total processing: ~200-300K cycles (~0.7-1.1 ms at 280 MHz)
  // - ADC chirp time: 64 µs (128 samples at 2 MSPS)
  // - 32 chirps: 2.048 ms → plenty of time for FFT processing
  //
  // VERIFICATION:
  // - adc_changes: Number of ADC samples captured (should be > 0)
  // - adc_value_before/after: Verify DMA is updating buffer
  //
  // INPUT DATA:
  // - adc_buffer[8192]: Raw 16-bit ADC samples (circular buffer)
  // - chirp_data[32][128]: Converted float32_t chirps from ADC data
  // - hamming_window_range[128]: Pre-computed 128-point Hamming window for Range FFT
  // - hamming_window_doppler[32]: Pre-computed 32-point Hamming window for Doppler FFT
  //
  // RANGE FFT RESULTS:
  // - fft_outputs[32][128]: 32 range FFT results organized by chirp (windowed)
  // - range_bins[64][64]: reorganized as 64 range bins x 32 complex samples
  //   * range_bins[0][0..63] = bin 0 from all 32 chirps (Re0,Im0, Re1,Im1, ..., Re31,Im31)
  //
  // DOPPLER FFT RESULTS:
  // - doppler_fft[64][64]: 64 Doppler FFT results (64 range bins x 32 velocity bins, windowed)
  //   * doppler_fft[range][0..63] = 32 complex velocity bins for this range bin
  //   * This is the final Range-Doppler map!
  //
  // PERFORMANCE METRICS:
  //
  // Range FFT Processing (32 chirps):
  // - total_range_window_cycles: total CPU cycles for 32 range window applications
  // - range_window_cycles: average cycles per range window (~100-300 expected)
  // - range_window_time_us: average range window time in microseconds
  // - total_range_fft_cycles: total CPU cycles for 32 range FFTs
  // - fft_cycles: average cycles per range FFT (~4-6K expected for 128-point)
  // - fft_time_us: average range FFT time in microseconds
  //
  // Doppler FFT Processing (64 range bins):
  // - total_doppler_window_cycles: total CPU cycles for 64 Doppler window applications
  // - doppler_window_cycles: average cycles per Doppler window (~50-100 expected)
  // - doppler_window_time_us: average Doppler window time in microseconds
  // - total_doppler_fft_cycles: total CPU cycles for 64 Doppler FFTs
  // - doppler_fft_cycles: average cycles per Doppler FFT (~2-3K expected)
  // - doppler_fft_time_us: average Doppler FFT time in microseconds
  //
  // RANGE-DOPPLER MAP PEAK DETECTION (5×5 Windowed Energy - HIGHLY OPTIMIZED):
  // - rd_map_peak_energy: Maximum windowed energy (sum of 25 cells in 5×5 window)
  // - rd_map_peak_range_bin: Range bin index (2-61) at center of peak window
  // - rd_map_peak_velocity_bin: Velocity bin index (2-29) at center of peak window
  // - rd_map_peak_range_m: Calculated range in meters (based on radar parameters)
  // - rd_map_peak_velocity_mps: Calculated velocity in m/s (based on radar parameters)
  // - peak_detection_cycles: CPU cycles for optimized search (~15-30K expected, was 550K before)
  // - peak_detection_time_us: Peak detection time in microseconds (~54-107 µs at 280 MHz)
  //
  // OPTIMIZATIONS APPLIED:
  // 1. Magnitude² with manual unrolling by 4: Better SIMD utilization
  // 2. Separable filtering (1×5 then 5×1): 8 ops/cell vs 24 (3× reduction)
  // 3. Sliding window in both passes: Only 2 ops per position after initialization
  // 4. Restrict pointers: Enables compiler auto-vectorization
  // 5. Cache-friendly access patterns: Sequential reads, minimal cache misses
  // - Expected speedup: 18-37× faster than naive approach (550K → 15-30K cycles)
  //
  // This represents the strongest radar target in both range and velocity dimensions.
  // The 5×5 window integrates energy over a neighborhood for better SNR and noise rejection.
  // Set a breakpoint here to inspect the peak location and energy.
  //
  // - fft_magnitude[64]: magnitude spectrum of first chirp (64 bins for 128-point real FFT)
  // - peak_frequency: depends on input signal

  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  {
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
    HAL_GPIO_TogglePin(USER_LED1_GPIO_Port, USER_LED1_Pin);
    HAL_Delay(1000);
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /*AXI clock gating */
  RCC->CKGAENR = 0xE003FFFF;

  /** Supply configuration update enable
  */
  HAL_PWREx_ConfigSupply(PWR_DIRECT_SMPS_SUPPLY);

  /** Configure the main internal regulator output voltage
  */
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE0);

  while(!__HAL_PWR_GET_FLAG(PWR_FLAG_VOSRDY)) {}

  /** Configure LSE Drive Capability
  */
  HAL_PWR_EnableBkUpAccess();
  __HAL_RCC_LSEDRIVE_CONFIG(RCC_LSEDRIVE_LOW);

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI|RCC_OSCILLATORTYPE_HSE
                              |RCC_OSCILLATORTYPE_LSE;
  RCC_OscInitStruct.HSEState = RCC_HSE_ON;
  RCC_OscInitStruct.LSEState = RCC_LSE_ON;
  RCC_OscInitStruct.HSIState = RCC_HSI_DIV1;
  RCC_OscInitStruct.HSICalibrationValue = 64;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSE;
  RCC_OscInitStruct.PLL.PLLM = 12;
  RCC_OscInitStruct.PLL.PLLN = 280;
  RCC_OscInitStruct.PLL.PLLP = 2;
  RCC_OscInitStruct.PLL.PLLQ = 3;
  RCC_OscInitStruct.PLL.PLLR = 4;
  RCC_OscInitStruct.PLL.PLLRGE = RCC_PLL1VCIRANGE_1;
  RCC_OscInitStruct.PLL.PLLVCOSEL = RCC_PLL1VCOWIDE;
  RCC_OscInitStruct.PLL.PLLFRACN = 0;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2
                              |RCC_CLOCKTYPE_D3PCLK1|RCC_CLOCKTYPE_D1PCLK1;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.SYSCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_HCLK_DIV2;
  RCC_ClkInitStruct.APB3CLKDivider = RCC_APB3_DIV2;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_APB1_DIV2;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_APB2_DIV2;
  RCC_ClkInitStruct.APB4CLKDivider = RCC_APB4_DIV2;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_3) != HAL_OK)
  {
    Error_Handler();
  }
  HAL_RCC_MCOConfig(RCC_MCO1, RCC_MCO1SOURCE_HSE, RCC_MCODIV_1);
}

/**
  * @brief Peripherals Common Clock Configuration
  * @retval None
  */
void PeriphCommonClock_Config(void)
{
  RCC_PeriphCLKInitTypeDef PeriphClkInitStruct = {0};

  /** Initializes the peripherals clock
  */
  PeriphClkInitStruct.PeriphClockSelection = RCC_PERIPHCLK_FMC|RCC_PERIPHCLK_OSPI
                              |RCC_PERIPHCLK_ADC|RCC_PERIPHCLK_CKPER;
  PeriphClkInitStruct.PLL2.PLL2M = 2;
  PeriphClkInitStruct.PLL2.PLL2N = 12;
  PeriphClkInitStruct.PLL2.PLL2P = 4;
  PeriphClkInitStruct.PLL2.PLL2Q = 4;
  PeriphClkInitStruct.PLL2.PLL2R = 4;
  PeriphClkInitStruct.PLL2.PLL2RGE = RCC_PLL2VCIRANGE_3;
  PeriphClkInitStruct.PLL2.PLL2VCOSEL = RCC_PLL2VCOWIDE;
  PeriphClkInitStruct.PLL2.PLL2FRACN = 0;
  PeriphClkInitStruct.FmcClockSelection = RCC_FMCCLKSOURCE_PLL2;
  PeriphClkInitStruct.OspiClockSelection = RCC_OSPICLKSOURCE_PLL2;
  PeriphClkInitStruct.CkperClockSelection = RCC_CLKPSOURCE_HSI;
  PeriphClkInitStruct.AdcClockSelection = RCC_ADCCLKSOURCE_PLL2;
  if (HAL_RCCEx_PeriphCLKConfig(&PeriphClkInitStruct) != HAL_OK)
  {
    Error_Handler();
  }
}

/**
  * @brief ADC1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_ADC1_Init(void)
{

  /* USER CODE BEGIN ADC1_Init 0 */

  /* USER CODE END ADC1_Init 0 */

  ADC_MultiModeTypeDef multimode = {0};
  ADC_ChannelConfTypeDef sConfig = {0};

  /* USER CODE BEGIN ADC1_Init 1 */

  /* USER CODE END ADC1_Init 1 */

  /** Common config
  */
  hadc1.Instance = ADC1;
  hadc1.Init.ClockPrescaler = ADC_CLOCK_ASYNC_DIV1;
  hadc1.Init.Resolution = ADC_RESOLUTION_16B;
  hadc1.Init.ScanConvMode = ADC_SCAN_DISABLE;
  hadc1.Init.EOCSelection = ADC_EOC_SINGLE_CONV;
  hadc1.Init.LowPowerAutoWait = DISABLE;
  hadc1.Init.ContinuousConvMode = ENABLE;
  hadc1.Init.NbrOfConversion = 1;
  hadc1.Init.DiscontinuousConvMode = DISABLE;
  hadc1.Init.ExternalTrigConv = ADC_SOFTWARE_START;
  hadc1.Init.ExternalTrigConvEdge = ADC_EXTERNALTRIGCONVEDGE_NONE;
  hadc1.Init.ConversionDataManagement = ADC_CONVERSIONDATA_DMA_CIRCULAR;
  hadc1.Init.Overrun = ADC_OVR_DATA_PRESERVED;
  hadc1.Init.LeftBitShift = ADC_LEFTBITSHIFT_NONE;
  hadc1.Init.OversamplingMode = DISABLE;
  if (HAL_ADC_Init(&hadc1) != HAL_OK)
  {
    Error_Handler();
  }

  /** Configure the ADC multi-mode
  */
  multimode.Mode = ADC_MODE_INDEPENDENT;
  if (HAL_ADCEx_MultiModeConfigChannel(&hadc1, &multimode) != HAL_OK)
  {
    Error_Handler();
  }

  /** Configure Regular Channel
  */
  sConfig.Channel = ADC_CHANNEL_1;
  sConfig.Rank = ADC_REGULAR_RANK_1;
  sConfig.SamplingTime = ADC_SAMPLETIME_1CYCLE_5;
  sConfig.SingleDiff = ADC_SINGLE_ENDED;
  sConfig.OffsetNumber = ADC_OFFSET_NONE;
  sConfig.Offset = 0;
  sConfig.OffsetSignedSaturation = DISABLE;
  if (HAL_ADC_ConfigChannel(&hadc1, &sConfig) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN ADC1_Init 2 */

  /* USER CODE END ADC1_Init 2 */

}

/**
  * @brief I2C4 Initialization Function
  * @param None
  * @retval None
  */
static void MX_I2C4_Init(void)
{

  /* USER CODE BEGIN I2C4_Init 0 */

  /* USER CODE END I2C4_Init 0 */

  /* USER CODE BEGIN I2C4_Init 1 */

  /* USER CODE END I2C4_Init 1 */
  hi2c4.Instance = I2C4;
  hi2c4.Init.Timing = 0xC010151E;
  hi2c4.Init.OwnAddress1 = 0;
  hi2c4.Init.AddressingMode = I2C_ADDRESSINGMODE_7BIT;
  hi2c4.Init.DualAddressMode = I2C_DUALADDRESS_DISABLE;
  hi2c4.Init.OwnAddress2 = 0;
  hi2c4.Init.OwnAddress2Masks = I2C_OA2_NOMASK;
  hi2c4.Init.GeneralCallMode = I2C_GENERALCALL_DISABLE;
  hi2c4.Init.NoStretchMode = I2C_NOSTRETCH_DISABLE;
  if (HAL_I2C_Init(&hi2c4) != HAL_OK)
  {
    Error_Handler();
  }

  /** Configure Analogue filter
  */
  if (HAL_I2CEx_ConfigAnalogFilter(&hi2c4, I2C_ANALOGFILTER_ENABLE) != HAL_OK)
  {
    Error_Handler();
  }

  /** Configure Digital filter
  */
  if (HAL_I2CEx_ConfigDigitalFilter(&hi2c4, 0) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN I2C4_Init 2 */

  /* USER CODE END I2C4_Init 2 */

}

/**
  * @brief I2S6 Initialization Function
  * @param None
  * @retval None
  */
static void MX_I2S6_Init(void)
{

  /* USER CODE BEGIN I2S6_Init 0 */

  /* USER CODE END I2S6_Init 0 */

  /* USER CODE BEGIN I2S6_Init 1 */

  /* USER CODE END I2S6_Init 1 */
  hi2s6.Instance = SPI6;
  hi2s6.Init.Mode = I2S_MODE_MASTER_FULLDUPLEX;
  hi2s6.Init.Standard = I2S_STANDARD_PHILIPS;
  hi2s6.Init.DataFormat = I2S_DATAFORMAT_16B;
  hi2s6.Init.MCLKOutput = I2S_MCLKOUTPUT_ENABLE;
  hi2s6.Init.AudioFreq = I2S_AUDIOFREQ_8K;
  hi2s6.Init.CPOL = I2S_CPOL_LOW;
  hi2s6.Init.FirstBit = I2S_FIRSTBIT_MSB;
  hi2s6.Init.WSInversion = I2S_WS_INVERSION_DISABLE;
  hi2s6.Init.Data24BitAlignment = I2S_DATA_24BIT_ALIGNMENT_RIGHT;
  hi2s6.Init.MasterKeepIOState = I2S_MASTER_KEEP_IO_STATE_DISABLE;
  if (HAL_I2S_Init(&hi2s6) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN I2S6_Init 2 */

  /* USER CODE END I2S6_Init 2 */

}

/**
  * @brief LTDC Initialization Function
  * @param None
  * @retval None
  */
static void MX_LTDC_Init(void)
{

  /* USER CODE BEGIN LTDC_Init 0 */

  /* USER CODE END LTDC_Init 0 */

  LTDC_LayerCfgTypeDef pLayerCfg = {0};
  LTDC_LayerCfgTypeDef pLayerCfg1 = {0};

  /* USER CODE BEGIN LTDC_Init 1 */

  /* USER CODE END LTDC_Init 1 */
  hltdc.Instance = LTDC;
  hltdc.Init.HSPolarity = LTDC_HSPOLARITY_AL;
  hltdc.Init.VSPolarity = LTDC_VSPOLARITY_AL;
  hltdc.Init.DEPolarity = LTDC_DEPOLARITY_AL;
  hltdc.Init.PCPolarity = LTDC_PCPOLARITY_IPC;
  hltdc.Init.HorizontalSync = 0;
  hltdc.Init.VerticalSync = 9;
  hltdc.Init.AccumulatedHBP = 43;
  hltdc.Init.AccumulatedVBP = 21;
  hltdc.Init.AccumulatedActiveW = 523;
  hltdc.Init.AccumulatedActiveH = 293;
  hltdc.Init.TotalWidth = 531;
  hltdc.Init.TotalHeigh = 297;
  hltdc.Init.Backcolor.Blue = 0;
  hltdc.Init.Backcolor.Green = 0;
  hltdc.Init.Backcolor.Red = 0;
  if (HAL_LTDC_Init(&hltdc) != HAL_OK)
  {
    Error_Handler();
  }
  pLayerCfg.WindowX0 = 0;
  pLayerCfg.WindowX1 = 0;
  pLayerCfg.WindowY0 = 0;
  pLayerCfg.WindowY1 = 0;
  pLayerCfg.PixelFormat = LTDC_PIXEL_FORMAT_ARGB8888;
  pLayerCfg.Alpha = 0;
  pLayerCfg.Alpha0 = 0;
  pLayerCfg.BlendingFactor1 = LTDC_BLENDING_FACTOR1_CA;
  pLayerCfg.BlendingFactor2 = LTDC_BLENDING_FACTOR2_CA;
  pLayerCfg.FBStartAdress = 0;
  pLayerCfg.ImageWidth = 0;
  pLayerCfg.ImageHeight = 0;
  pLayerCfg.Backcolor.Blue = 0;
  pLayerCfg.Backcolor.Green = 0;
  pLayerCfg.Backcolor.Red = 0;
  if (HAL_LTDC_ConfigLayer(&hltdc, &pLayerCfg, 0) != HAL_OK)
  {
    Error_Handler();
  }
  pLayerCfg1.WindowX0 = 0;
  pLayerCfg1.WindowX1 = 0;
  pLayerCfg1.WindowY0 = 0;
  pLayerCfg1.WindowY1 = 0;
  pLayerCfg1.PixelFormat = LTDC_PIXEL_FORMAT_ARGB8888;
  pLayerCfg1.Alpha = 0;
  pLayerCfg1.Alpha0 = 0;
  pLayerCfg1.BlendingFactor1 = LTDC_BLENDING_FACTOR1_CA;
  pLayerCfg1.BlendingFactor2 = LTDC_BLENDING_FACTOR2_CA;
  pLayerCfg1.FBStartAdress = 0;
  pLayerCfg1.ImageWidth = 0;
  pLayerCfg1.ImageHeight = 0;
  pLayerCfg1.Backcolor.Blue = 0;
  pLayerCfg1.Backcolor.Green = 0;
  pLayerCfg1.Backcolor.Red = 0;
  if (HAL_LTDC_ConfigLayer(&hltdc, &pLayerCfg1, 1) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN LTDC_Init 2 */

  /* USER CODE END LTDC_Init 2 */

}

/**
  * @brief OCTOSPI1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_OCTOSPI1_Init(void)
{

  /* USER CODE BEGIN OCTOSPI1_Init 0 */

  /* USER CODE END OCTOSPI1_Init 0 */

  OSPIM_CfgTypeDef sOspiManagerCfg = {0};

  /* USER CODE BEGIN OCTOSPI1_Init 1 */

  /* USER CODE END OCTOSPI1_Init 1 */
  /* OCTOSPI1 parameter configuration*/
  hospi1.Instance = OCTOSPI1;
  hospi1.Init.FifoThreshold = 1;
  hospi1.Init.DualQuad = HAL_OSPI_DUALQUAD_DISABLE;
  hospi1.Init.MemoryType = HAL_OSPI_MEMTYPE_MACRONIX;
  hospi1.Init.DeviceSize = 32;
  hospi1.Init.ChipSelectHighTime = 1;
  hospi1.Init.FreeRunningClock = HAL_OSPI_FREERUNCLK_DISABLE;
  hospi1.Init.ClockMode = HAL_OSPI_CLOCK_MODE_0;
  hospi1.Init.WrapSize = HAL_OSPI_WRAP_NOT_SUPPORTED;
  hospi1.Init.ClockPrescaler = 1;
  hospi1.Init.SampleShifting = HAL_OSPI_SAMPLE_SHIFTING_NONE;
  hospi1.Init.DelayHoldQuarterCycle = HAL_OSPI_DHQC_DISABLE;
  hospi1.Init.ChipSelectBoundary = 0;
  hospi1.Init.DelayBlockBypass = HAL_OSPI_DELAY_BLOCK_BYPASSED;
  hospi1.Init.MaxTran = 0;
  hospi1.Init.Refresh = 0;
  if (HAL_OSPI_Init(&hospi1) != HAL_OK)
  {
    Error_Handler();
  }
  sOspiManagerCfg.ClkPort = 1;
  sOspiManagerCfg.DQSPort = 1;
  sOspiManagerCfg.NCSPort = 1;
  sOspiManagerCfg.IOLowPort = HAL_OSPIM_IOPORT_1_HIGH;
  if (HAL_OSPIM_Config(&hospi1, &sOspiManagerCfg, HAL_OSPI_TIMEOUT_DEFAULT_VALUE) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN OCTOSPI1_Init 2 */

  /* USER CODE END OCTOSPI1_Init 2 */

}

/**
  * @brief RTC Initialization Function
  * @param None
  * @retval None
  */
static void MX_RTC_Init(void)
{

  /* USER CODE BEGIN RTC_Init 0 */

  /* USER CODE END RTC_Init 0 */

  /* USER CODE BEGIN RTC_Init 1 */

  /* USER CODE END RTC_Init 1 */

  /** Initialize RTC Only
  */
  hrtc.Instance = RTC;
  hrtc.Init.HourFormat = RTC_HOURFORMAT_24;
  hrtc.Init.AsynchPrediv = 127;
  hrtc.Init.SynchPrediv = 255;
  hrtc.Init.OutPut = RTC_OUTPUT_DISABLE;
  hrtc.Init.OutPutPolarity = RTC_OUTPUT_POLARITY_HIGH;
  hrtc.Init.OutPutType = RTC_OUTPUT_TYPE_OPENDRAIN;
  hrtc.Init.OutPutRemap = RTC_OUTPUT_REMAP_NONE;
  if (HAL_RTC_Init(&hrtc) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN RTC_Init 2 */

  /* USER CODE END RTC_Init 2 */

}

/**
  * @brief SPI2 Initialization Function
  * @param None
  * @retval None
  */
static void MX_SPI2_Init(void)
{

  /* USER CODE BEGIN SPI2_Init 0 */

  /* USER CODE END SPI2_Init 0 */

  /* USER CODE BEGIN SPI2_Init 1 */

  /* USER CODE END SPI2_Init 1 */
  /* SPI2 parameter configuration*/
  hspi2.Instance = SPI2;
  hspi2.Init.Mode = SPI_MODE_MASTER;
  hspi2.Init.Direction = SPI_DIRECTION_2LINES;
  hspi2.Init.DataSize = SPI_DATASIZE_4BIT;
  hspi2.Init.CLKPolarity = SPI_POLARITY_LOW;
  hspi2.Init.CLKPhase = SPI_PHASE_1EDGE;
  hspi2.Init.NSS = SPI_NSS_SOFT;
  hspi2.Init.BaudRatePrescaler = SPI_BAUDRATEPRESCALER_2;
  hspi2.Init.FirstBit = SPI_FIRSTBIT_MSB;
  hspi2.Init.TIMode = SPI_TIMODE_DISABLE;
  hspi2.Init.CRCCalculation = SPI_CRCCALCULATION_DISABLE;
  hspi2.Init.CRCPolynomial = 0x0;
  hspi2.Init.NSSPMode = SPI_NSS_PULSE_ENABLE;
  hspi2.Init.NSSPolarity = SPI_NSS_POLARITY_LOW;
  hspi2.Init.FifoThreshold = SPI_FIFO_THRESHOLD_01DATA;
  hspi2.Init.TxCRCInitializationPattern = SPI_CRC_INITIALIZATION_ALL_ZERO_PATTERN;
  hspi2.Init.RxCRCInitializationPattern = SPI_CRC_INITIALIZATION_ALL_ZERO_PATTERN;
  hspi2.Init.MasterSSIdleness = SPI_MASTER_SS_IDLENESS_00CYCLE;
  hspi2.Init.MasterInterDataIdleness = SPI_MASTER_INTERDATA_IDLENESS_00CYCLE;
  hspi2.Init.MasterReceiverAutoSusp = SPI_MASTER_RX_AUTOSUSP_DISABLE;
  hspi2.Init.MasterKeepIOState = SPI_MASTER_KEEP_IO_STATE_DISABLE;
  hspi2.Init.IOSwap = SPI_IO_SWAP_DISABLE;
  if (HAL_SPI_Init(&hspi2) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN SPI2_Init 2 */

  /* USER CODE END SPI2_Init 2 */

}

/**
  * @brief USART1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_USART1_UART_Init(void)
{

  /* USER CODE BEGIN USART1_Init 0 */

  /* USER CODE END USART1_Init 0 */

  /* USER CODE BEGIN USART1_Init 1 */

  /* USER CODE END USART1_Init 1 */
  huart1.Instance = USART1;
  huart1.Init.BaudRate = 115200;
  huart1.Init.WordLength = UART_WORDLENGTH_8B;
  huart1.Init.StopBits = UART_STOPBITS_1;
  huart1.Init.Parity = UART_PARITY_NONE;
  huart1.Init.Mode = UART_MODE_TX_RX;
  huart1.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  huart1.Init.OverSampling = UART_OVERSAMPLING_16;
  huart1.Init.OneBitSampling = UART_ONE_BIT_SAMPLE_DISABLE;
  huart1.Init.ClockPrescaler = UART_PRESCALER_DIV1;
  huart1.AdvancedInit.AdvFeatureInit = UART_ADVFEATURE_NO_INIT;
  if (HAL_UART_Init(&huart1) != HAL_OK)
  {
    Error_Handler();
  }
  if (HAL_UARTEx_SetTxFifoThreshold(&huart1, UART_TXFIFO_THRESHOLD_1_8) != HAL_OK)
  {
    Error_Handler();
  }
  if (HAL_UARTEx_SetRxFifoThreshold(&huart1, UART_RXFIFO_THRESHOLD_1_8) != HAL_OK)
  {
    Error_Handler();
  }
  if (HAL_UARTEx_DisableFifoMode(&huart1) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN USART1_Init 2 */

  /* USER CODE END USART1_Init 2 */

}

/**
  * Enable DMA controller clock
  */
static void MX_DMA_Init(void)
{

  /* DMA controller clock enable */
  __HAL_RCC_DMA1_CLK_ENABLE();

  /* DMA interrupt init */
  /* DMA1_Stream0_IRQn interrupt configuration */
  HAL_NVIC_SetPriority(DMA1_Stream0_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(DMA1_Stream0_IRQn);

}

/* FMC initialization function */
static void MX_FMC_Init(void)
{

  /* USER CODE BEGIN FMC_Init 0 */

  /* USER CODE END FMC_Init 0 */

  FMC_SDRAM_TimingTypeDef SdramTiming = {0};

  /* USER CODE BEGIN FMC_Init 1 */

  /* USER CODE END FMC_Init 1 */

  /** Perform the SDRAM1 memory initialization sequence
  */
  hsdram1.Instance = FMC_SDRAM_DEVICE;
  /* hsdram1.Init */
  hsdram1.Init.SDBank = FMC_SDRAM_BANK2;
  hsdram1.Init.ColumnBitsNumber = FMC_SDRAM_COLUMN_BITS_NUM_8;
  hsdram1.Init.RowBitsNumber = FMC_SDRAM_ROW_BITS_NUM_12;
  hsdram1.Init.MemoryDataWidth = FMC_SDRAM_MEM_BUS_WIDTH_16;
  hsdram1.Init.InternalBankNumber = FMC_SDRAM_INTERN_BANKS_NUM_4;
  hsdram1.Init.CASLatency = FMC_SDRAM_CAS_LATENCY_1;
  hsdram1.Init.WriteProtection = FMC_SDRAM_WRITE_PROTECTION_DISABLE;
  hsdram1.Init.SDClockPeriod = FMC_SDRAM_CLOCK_DISABLE;
  hsdram1.Init.ReadBurst = FMC_SDRAM_RBURST_DISABLE;
  hsdram1.Init.ReadPipeDelay = FMC_SDRAM_RPIPE_DELAY_0;
  /* SdramTiming */
  SdramTiming.LoadToActiveDelay = 16;
  SdramTiming.ExitSelfRefreshDelay = 16;
  SdramTiming.SelfRefreshTime = 16;
  SdramTiming.RowCycleDelay = 16;
  SdramTiming.WriteRecoveryTime = 16;
  SdramTiming.RPDelay = 16;
  SdramTiming.RCDDelay = 16;

  if (HAL_SDRAM_Init(&hsdram1, &SdramTiming) != HAL_OK)
  {
    Error_Handler( );
  }

  /* USER CODE BEGIN FMC_Init 2 */

  /* USER CODE END FMC_Init 2 */
}

/**
  * @brief GPIO Initialization Function
  * @param None
  * @retval None
  */
static void MX_GPIO_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct = {0};
  /* USER CODE BEGIN MX_GPIO_Init_1 */

  /* USER CODE END MX_GPIO_Init_1 */

  /* GPIO Ports Clock Enable */
  __HAL_RCC_GPIOI_CLK_ENABLE();
  __HAL_RCC_GPIOG_CLK_ENABLE();
  __HAL_RCC_GPIOK_CLK_ENABLE();
  __HAL_RCC_GPIOD_CLK_ENABLE();
  __HAL_RCC_GPIOC_CLK_ENABLE();
  __HAL_RCC_GPIOE_CLK_ENABLE();
  __HAL_RCC_GPIOJ_CLK_ENABLE();
  __HAL_RCC_GPIOA_CLK_ENABLE();
  __HAL_RCC_GPIOB_CLK_ENABLE();
  __HAL_RCC_GPIOF_CLK_ENABLE();
  __HAL_RCC_GPIOH_CLK_ENABLE();

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOI, WIFI_BOOT_Pin|WIFI_WKUP_Pin|WIFI_RST_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOG, USER_LED1_Pin|USER_LED2_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(SPI2_NSS_GPIO_Port, SPI2_NSS_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(AUDIO_NRST_GPIO_Port, AUDIO_NRST_Pin, GPIO_PIN_SET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOA, LCD_BL_CTRL_Pin|LCD_ON_OFF_Pin, GPIO_PIN_SET);

  /*Configure GPIO pins : WIFI_GPIO_Pin WIFI_DATRDY_Pin */
  GPIO_InitStruct.Pin = WIFI_GPIO_Pin|WIFI_DATRDY_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(GPIOI, &GPIO_InitStruct);

  /*Configure GPIO pins : SDIO1_D2_Pin SDIO1_CK_Pin SDIO1_D3_Pin SDIO1_D1_Pin
                           SDIO1_D0_Pin */
  GPIO_InitStruct.Pin = SDIO1_D2_Pin|SDIO1_CK_Pin|SDIO1_D3_Pin|SDIO1_D1_Pin
                          |SDIO1_D0_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_VERY_HIGH;
  GPIO_InitStruct.Alternate = GPIO_AF12_SDMMC1;
  HAL_GPIO_Init(GPIOC, &GPIO_InitStruct);

  /*Configure GPIO pins : WIFI_BOOT_Pin WIFI_WKUP_Pin WIFI_RST_Pin */
  GPIO_InitStruct.Pin = WIFI_BOOT_Pin|WIFI_WKUP_Pin|WIFI_RST_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOI, &GPIO_InitStruct);

  /*Configure GPIO pins : USER_LED1_Pin AUDIO_NRST_Pin USER_LED2_Pin */
  GPIO_InitStruct.Pin = USER_LED1_Pin|AUDIO_NRST_Pin|USER_LED2_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOG, &GPIO_InitStruct);

  /*Configure GPIO pin : SDIO1_CMD_Pin */
  GPIO_InitStruct.Pin = SDIO1_CMD_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_VERY_HIGH;
  GPIO_InitStruct.Alternate = GPIO_AF12_SDMMC1;
  HAL_GPIO_Init(SDIO1_CMD_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : uSD_Detect_Pin */
  GPIO_InitStruct.Pin = uSD_Detect_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStruct.Pull = GPIO_PULLUP;
  HAL_GPIO_Init(uSD_Detect_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pins : SPI2_NSS_Pin LCD_BL_CTRL_Pin LCD_ON_OFF_Pin */
  GPIO_InitStruct.Pin = SPI2_NSS_Pin|LCD_BL_CTRL_Pin|LCD_ON_OFF_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOA, &GPIO_InitStruct);

  /*Configure GPIO pin : WAKEUP_Pin */
  GPIO_InitStruct.Pin = WAKEUP_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(WAKEUP_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : MCO_Pin */
  GPIO_InitStruct.Pin = MCO_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  GPIO_InitStruct.Alternate = GPIO_AF0_MCO;
  HAL_GPIO_Init(MCO_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : LCD_INT_Pin */
  GPIO_InitStruct.Pin = LCD_INT_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(LCD_INT_GPIO_Port, &GPIO_InitStruct);

  /*AnalogSwitch Config */
  HAL_SYSCFG_AnalogSwitchConfig(SYSCFG_SWITCH_PA1, SYSCFG_SWITCH_PA1_CLOSE);

  /* USER CODE BEGIN MX_GPIO_Init_2 */

  /* USER CODE END MX_GPIO_Init_2 */
}

/* USER CODE BEGIN 4 */

/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}
#ifdef USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
