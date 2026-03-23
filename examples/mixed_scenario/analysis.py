import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.signal import convolve, get_window, stft, decimate
import os

# --- Scenario Constants ---
C = 299792458.0
FS = 10.0e6  # ADC Sample Rate (10 MHz)

# Target Kinematics
TGT_START_POS = np.array([500.0, 100.0, 100.0])
TGT_VELOCITY = np.array([0.0, -200.0, 0.0])  # Moves Y: 100 to -100 in 1 sec (Mach ~0.58)

# Pulsed Radar Parameters
PULSED_FILE = "PulsedRadar_results.h5"
PULSE_REF_FILE = "pulse.h5"
PULSED_POS = np.array([0.0, 0.0, 0.0])
FC_PULSED = 10.0e9
LAMBDA_PULSED = C / FC_PULSED
PRF = 10000.0
WINDOW_SKIP = 1e-6

# CW Radar Parameters
CW_FILE = "CWRadar_results.h5"
CW_POS = np.array([1000.0, 0.0, 0.0])
FC_CW = 5.0e9
LAMBDA_CW = C / FC_CW
PT_CW = 1000.0
RCS = 1.0


# ==========================================
# HELPER FUNCTIONS
# ==========================================

def get_target_state(t):
    """Returns target position and velocity at time t."""
    return TGT_START_POS + (TGT_VELOCITY * t), TGT_VELOCITY


def calc_monostatic_theory(t, radar_pos, wavelength):
    """Calculates theoretical monostatic Range and Doppler."""
    tgt_pos, tgt_vel = get_target_state(t)

    vec_radar_tgt = tgt_pos - radar_pos
    rng = np.linalg.norm(vec_radar_tgt)

    # Radial velocity (positive = moving away)
    u_vec = vec_radar_tgt / rng
    v_rad = np.dot(tgt_vel, u_vec)

    # Doppler (FERS Convention: Approaching = Positive Shift)
    doppler = -2.0 * v_rad / wavelength

    return rng, doppler


# ==========================================
# PULSED RADAR ANALYSIS
# ==========================================

def analyze_pulsed():
    print("\n" + "=" * 50)
    print("--- STARTING PULSED RADAR ANALYSIS ---")
    print("=" * 50)

    if not os.path.exists(PULSED_FILE) or not os.path.exists(PULSE_REF_FILE):
        print(f"ERROR: Missing {PULSED_FILE} or {PULSE_REF_FILE}")
        return

    # 1. Load Reference Pulse
    with h5py.File(PULSE_REF_FILE, 'r') as f:
        ref_pulse = f['I/value'][:] + 1j * f['Q/value'][:]

    # 2. Load FERS Data
    pulses = []
    with h5py.File(PULSED_FILE, 'r') as f:
        chunk_names = sorted([k for k in f.keys() if k.endswith('_I')])
        for key_i in chunk_names:
            key_q = key_i.replace('_I', '_Q')
            fullscale = f[key_i].attrs['fullscale']
            complex_pulse = (f[key_i][:] + 1j * f[key_q][:]) * fullscale
            pulses.append(complex_pulse)

    raw_data = np.array(pulses)

    # Process a Coherent Processing Interval (CPI)
    cpi_len = min(256, raw_data.shape[0])
    cpi_data = raw_data[:cpi_len, :]
    print(f"Processing CPI: {cpi_len} pulses")

    # 3. Pulse Compression & Doppler FFT
    matched_filter = np.conj(ref_pulse[::-1])
    compressed = np.array([convolve(p, matched_filter, mode='full') for p in cpi_data])

    win = get_window('hamming', cpi_len)
    windowed_matrix = compressed * win[:, np.newaxis]
    rd_map = np.fft.fftshift(np.fft.fft(windowed_matrix, axis=0), axes=0)
    rd_map_mag = np.abs(rd_map)

    # 4. Axis Construction
    filter_lag = len(ref_pulse) - 1
    corrected_indices = np.arange(rd_map.shape[1]) - filter_lag
    fast_time_axis = WINDOW_SKIP + (corrected_indices / FS)
    range_axis = fast_time_axis * C / 2.0  # Monostatic: divide by 2

    doppler_axis = np.fft.fftshift(np.fft.fftfreq(cpi_len, d=1 / PRF))

    # 5. Find Peak
    peak_idx_dop, peak_idx_rng = np.unravel_index(np.argmax(rd_map_mag), rd_map_mag.shape)
    measured_range = range_axis[peak_idx_rng]
    measured_doppler = doppler_axis[peak_idx_dop]

    # 6. Theoretical Verification
    t_mid = (cpi_len / 2) / PRF
    theo_range, theo_doppler = calc_monostatic_theory(t_mid, PULSED_POS, LAMBDA_PULSED)

    # Target is moving at Mach 0.58 (200 m/s). Max radial velocity is ~38.5 m/s.
    # Max Doppler is ~2.56 kHz. Since PRF is 10 kHz, Doppler is NOT aliased.
    theo_doppler_aliased = (theo_doppler + PRF / 2) % PRF - PRF / 2

    print(f"{'Metric':<20} | {'Theoretical':<15} | {'Measured':<15} | {'Error':<15}")
    print("-" * 70)
    print(
        f"{'Range (m)':<20} | {theo_range:<15.4f} | {measured_range:<15.4f} | {abs(measured_range - theo_range):<15.4f}")
    print(
        f"{'Doppler (Hz)':<20} | {theo_doppler_aliased:<15.4f} | {measured_doppler:<15.4f} | {abs(measured_doppler - theo_doppler_aliased):<15.4f}")
    print("-" * 70)

    # 7. Plotting
    plt.figure(figsize=(12, 5))

    plt.subplot(1, 2, 1)
    rd_db = 20 * np.log10(rd_map_mag + 1e-9)
    rd_db -= np.max(rd_db)
    extent = [range_axis[0], range_axis[-1], doppler_axis[0], doppler_axis[-1]]
    plt.imshow(rd_db, aspect='auto', extent=extent, origin='lower', cmap='jet', vmin=-60)
    plt.scatter([theo_range], [theo_doppler_aliased], color='white', marker='x', s=100, label='Theory')
    plt.title('Pulsed Radar: Range-Doppler Map')
    plt.xlabel('Range (m)')
    plt.ylabel('Doppler (Hz)')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(range_axis, rd_map_mag[peak_idx_dop, :], label='Measured Profile')
    plt.axvline(theo_range, color='r', linestyle='--', label='Theory Range')
    plt.xlim(max(0, theo_range - 300), theo_range + 300)
    plt.title('Range Profile (Cut at Peak Doppler)')
    plt.xlabel('Range (m)')
    plt.ylabel('Magnitude')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.savefig("pulsed_analysis.png", dpi=300)
    print(">> Saved pulsed_analysis.png")


# ==========================================
# CW RADAR ANALYSIS
# ==========================================

def analyze_cw():
    print("\n" + "=" * 50)
    print("--- STARTING CW RADAR ANALYSIS ---")
    print("=" * 50)

    if not os.path.exists(CW_FILE):
        print(f"ERROR: Missing {CW_FILE}")
        return

    with h5py.File(CW_FILE, 'r') as f:
        fs = f.attrs['sampling_rate']
        start_time = f.attrs['start_time']
        fullscale = f.attrs.get('fullscale', 1.0)
        iq_data = (f['I_data'][:] + 1j * f['Q_data'][:]) * fullscale

    num_samples = len(iq_data)
    time_vector = np.linspace(start_time, start_time + (num_samples / fs), num_samples, endpoint=False)

    # 1. Theoretical Model
    theo_doppler = np.zeros(num_samples)
    theo_power = np.zeros(num_samples)
    for i, t in enumerate(time_vector):
        rng, dop = calc_monostatic_theory(t, CW_POS, LAMBDA_CW)
        theo_doppler[i] = dop
        theo_power[i] = (PT_CW * (LAMBDA_CW ** 2) * RCS) / (((4 * np.pi) ** 3) * (rng ** 4))

    # 2. Decimation & STFT Processing
    # Processing 10M samples for CW is slow and yields poor frequency resolution.
    # We decimate by 100 to get a 100 kHz sample rate, allowing finer frequency bins.
    decimation_factor = 100
    fs_dec = fs / decimation_factor

    print(f"Decimating CW data from {fs / 1e6:.1f} MHz to {fs_dec / 1e3:.1f} kHz for STFT...")
    iq_dec = decimate(iq_data, decimation_factor)

    # STFT on decimated data (4096 bins @ 100kHz = ~24.4 Hz resolution)
    nperseg = 4096
    f_stft, t_stft, Zxx = stft(iq_dec, fs_dec, nperseg=nperseg, noverlap=nperseg // 2, return_onesided=False)
    Zxx = np.fft.fftshift(Zxx, axes=0)
    f_stft = np.fft.fftshift(f_stft)

    # Adjust STFT time vector to match simulation start time
    t_stft += start_time

    peak_indices = np.argmax(np.abs(Zxx), axis=0)
    measured_doppler_stft = f_stft[peak_indices]
    theo_doppler_interp = np.interp(t_stft, time_vector, theo_doppler)
    doppler_mae = np.mean(np.abs(measured_doppler_stft - theo_doppler_interp))

    # Power Analysis (using original rate for accuracy)
    measured_power = np.abs(iq_data) ** 2
    power_rmse_db = np.sqrt(np.mean((10 * np.log10(measured_power + 1e-30) - 10 * np.log10(theo_power + 1e-30)) ** 2))

    print(f"CW Doppler MAE (STFT Peak): {doppler_mae:.2f} Hz")
    print(f"CW Power RMSE:              {power_rmse_db:.2f} dB")
    print(f"Doppler Sweep Range:        {np.min(theo_doppler):.1f} Hz to {np.max(theo_doppler):.1f} Hz")

    # 3. Plotting
    # Downsample time-domain power for plotting to avoid GUI lag
    plot_step = max(1, num_samples // 2000)
    time_plot = time_vector[::plot_step]
    power_plot = measured_power[::plot_step]
    theo_power_plot = theo_power[::plot_step]

    fig = plt.figure(figsize=(12, 10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1.5])

    # Plot 1: Power
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(time_plot, 10 * np.log10(power_plot + 1e-30), 'b-', alpha=0.7, label='Measured Power (Downsampled)')
    ax1.plot(time_plot, 10 * np.log10(theo_power_plot + 1e-30), 'r--', linewidth=2, label='Theoretical Power')
    ax1.set_title('CW Radar: Received Power')
    ax1.set_ylabel('Power (dBW)')
    ax1.set_xlabel('Time (s)')
    ax1.legend()
    ax1.grid(True)

    # Plot 2: STFT Spectrogram
    ax2 = fig.add_subplot(gs[1])
    Zxx_db = 20 * np.log10(np.abs(Zxx) + 1e-12)
    Zxx_db -= np.max(Zxx_db)

    # Zoom in on the relevant Doppler band
    zoom = 2000
    mask = (f_stft >= -zoom) & (f_stft <= zoom)

    mesh = ax2.pcolormesh(t_stft, f_stft[mask], Zxx_db[mask, :], shading='auto', cmap='inferno', vmin=-50, vmax=0)
    fig.colorbar(mesh, ax=ax2, label='Normalized Power (dB)')
    ax2.plot(time_vector, theo_doppler, 'c--', linewidth=2, label='Theoretical Doppler Track')
    ax2.set_title('CW Radar: High-Resolution STFT Spectrogram')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Frequency (Hz)')
    ax2.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig("cw_analysis.png", dpi=300)
    print(">> Saved cw_analysis.png")


if __name__ == "__main__":
    analyze_pulsed()
    analyze_cw()
    print("\nAnalysis Complete.")
