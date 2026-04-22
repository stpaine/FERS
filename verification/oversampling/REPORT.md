# Oversampling Audit Report: `libfers`

## 1. Executive Summary
An extensive audit of the `libfers` oversampling implementation was conducted across a matrix of oversampling ratios (1x to 64x) and filter lengths (8 to 64). 

The audit reveals that the core logic, clocking, and noise generation of the system are mathematically sound and bug-free. However, there is a **fundamental architectural limitation in the DSP filter design**. The system utilizes a single-stage Blackman-windowed sinc FIR filter. At high oversampling ratios, if the filter length is not scaled proportionally, the filter window truncates the sinc wave prematurely. This results in severe signal attenuation, phase distortion, and fractional delay errors.

Crucially, because `libfers` in production is **hardcoded to a filter length of 33**, the system is only mathematically stable for oversampling ratios up to **8x**. Ratios of 16x and above will result in silent data corruption.

---

## 2. Audit Results Overview
The audit executed 128 test cases across four suites, yielding an overall pass rate of **64.84%**.

| Suite | Pass Rate | Observation |
| :--- | :--- | :--- |
| **Clocking** | 100% | PRF quantization, skip logic, and CW buffer sizing are perfectly implemented. |
| **Pipeline** | 100% | Thermal noise generation (Boltzmann equations) and ADC quantization are flawless. |
| **Resampler** | 37.5% | Round-trip fidelity (upsample $\rightarrow$ downsample) fails predictably when the filter length is too short for the requested ratio. |
| **Render** | 21.8% | Highly sensitive to phase/amplitude ripple. Fails strict tolerances (`< 1e-3`) whenever the underlying FIR filter is starved for taps. |

---

## 3. Core Limitation: The "Practical Envelope"
To oversample a signal by a ratio of $N$ without losing energy, the sum of the FIR filter coefficients must equal $N$. As $N$ increases, the sinc wave stretches in the time domain. 

The audit establishes a **Practical Envelope** for this single-stage filter design:
> **`filter_length >= 4 * oversample_ratio`**

If this envelope is breached, the filter window cuts off the sinc wave before it accumulates enough area. 
* **Extreme Example:** At `ratio=64` and `filter_length=8`, the FIR sum is only `6.7` (instead of 64.0). The resulting signal retains only **~1%** of its original energy (`estimated_roundtrip_gain = 0.0109`).

---

## 4. Production Impact: The Hardcoded `filter_length = 33`
In production, `libfers` fixes the filter length at 33. Applying the Practical Envelope ($33 / 4 = 8.25$), the **maximum safe oversampling ratio is 8x**. 

Here is the exact impact on production users based on their configured ratio:

### The Safe Zone (Ratios 1x, 2x, 4x)
* **Status:** Perfect fidelity.
* **Impact:** The filter length of 33 is more than sufficient. Signal gain is preserved at 100%, and fractional delays are calculated with near-perfect precision.

### The Boundary (Ratio 8x)
* **Status:** Acceptable for production, but mathematically degraded.
* **Impact:** The system sits right on the edge of the envelope ($8 \times 4 = 32$). Round-trip gain is preserved (`1.0006`), but the system fails the audit's ultra-strict `< 1e-5` tolerance for constant rendering. Users will experience microscopic phase/amplitude ripples, though these are likely imperceptible in standard radar simulations.

### The Danger Zone (Ratios 16x, 24x, 32x, 64x)
* **Status:** **Silent Failure / Data Corruption.**
* **Impact:** Users configuring high ratios to achieve "higher fidelity" will actually destroy their simulation data. Because the filter length does not scale, the signal is massively attenuated.
    * **16x Ratio:** ~3.4% energy loss (`gain = 0.966`)
    * **24x Ratio:** ~25.0% energy loss (`gain = 0.751`)
    * **32x Ratio:** ~46.5% energy loss (`gain = 0.535`)
    * **64x Ratio:** ~82.9% energy loss (`gain = 0.171`)

---

## 5. Design Choices for Rectification
To prevent silent failures and improve the robustness of `libfers`, the following design choices are recommended, ordered from immediate mitigations to long-term architectural improvements.

### Option 1: Hard Cap & Fail-Fast (Low Effort)
If the filter length must remain hardcoded to 33 for performance or legacy reasons, the system must protect the user from invalid configurations.
* **Action:** Add a validation check during initialization.
* **Implementation:**
  ```cpp
  if (params::oversampleRatio() > 8) {
      throw std::invalid_argument(
          "Oversampling ratios > 8 are not supported with the current fixed filter length of 33. "
          "Signal attenuation will occur."
      );
  }
  ```

### Option 2: Dynamic Filter Length (Medium Effort)
Decouple the filter length from the hardcoded `33` and calculate it dynamically based on the requested ratio.
* **Action:** Update `params::renderFilterLength()` to scale with the ratio.
* **Implementation:**
  ```cpp
  unsigned calculateSafeFilterLength(unsigned ratio) {
      // Maintain legacy minimum of 33, but scale up for higher ratios
      return std::max(33u, 4 * ratio); 
  }
  ```
* **Trade-off:** A single-stage FIR filter's computational cost scales linearly with filter length. Oversampling by 64x with a filter length of 256 will cause a massive CPU performance hit during rendering.

### Option 3: Advisory Warnings (Minimum Viable Mitigation)
If backward compatibility dictates that high ratios *must* be allowed to run despite the attenuation, port the audit harness's warning logic into production.
* **Action:** Emit a `LOG(WARNING)` when the envelope is breached.
* **Implementation:**
  ```cpp
  if (filter_length < 4 * ratio) {
      LOG(WARNING, "Oversampling FIR under-spec'd for ratio=", ratio, 
                   ". Severe signal attenuation may occur.");
  }
  ```

### Option 4: Polyphase Interpolation (High Effort)
If `libfers` requires high-fidelity oversampling at ratios of 16x, 32x, or 64x, the current single-stage FIR design is mathematically and computationally inappropriate.
* **Action:** Replace the single-stage upsampler with a multi-stage polyphase interpolator.
* **Implementation:** To achieve 64x oversampling, cascade three 4x polyphase filters ($4 \times 4 \times 4 = 64$). 
* **Benefit:** This allows for short, highly efficient filter lengths at each stage, preserving 100% of the signal energy without the exponential CPU cost of a massive single-stage FIR.
