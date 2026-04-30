#!/usr/bin/env python3
import argparse
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = 299_792_458.0
RESULT_FILE = "BistaticFmcwRx_results.h5"
PLOT_FILE = "bistatic_external_dechirp_analysis.png"

FC = 5.0e9
CHIRP_BANDWIDTH = 1.0e6
CHIRP_DURATION = 1.0e-3
CHIRP_PERIOD = 1.0e-3
CHIRP_COUNT = 80
START_FREQUENCY_OFFSET = -5.0e5
CHIRP_RATE = CHIRP_BANDWIDTH / CHIRP_DURATION

TX_POS = np.array([-100.0, 0.0, 0.0])
RX_POS = np.array([100.0, 0.0, 0.0])
TARGET_START = np.array([0.0, 300.0, 0.0])
TARGET_VELOCITY = np.array([0.0, -10.0, 0.0])


def attr_text(attrs, name):
    value = attrs[name]
    return value.decode("utf-8") if isinstance(value, bytes) else value


def load_raw_iq(path):
    if not path.exists():
        raise FileNotFoundError(f"missing simulation output: {path}")

    with h5py.File(path, "r") as h5:
        fs = float(h5.attrs["sampling_rate"])
        start_time = float(h5.attrs["start_time"])
        fullscale = float(h5.attrs["fullscale"])
        metadata = json.loads(attr_text(h5.attrs, "fers_metadata_json"))
        iq = (h5["I_data"][:] + 1j * h5["Q_data"][:]) * fullscale

        assert attr_text(h5.attrs, "data_mode") == "fmcw"
        assert attr_text(h5.attrs, "fmcw_dechirp_mode") == "none"
        assert attr_text(h5.attrs, "fmcw_dechirp_reference_source") == "none"
        assert int(h5.attrs["fmcw_source_count"]) == 1
        assert int(h5.attrs["fmcw_chirp_count"]) == CHIRP_COUNT
        assert np.isclose(float(h5.attrs["fmcw_chirp_bandwidth"]), CHIRP_BANDWIDTH)
        assert np.isclose(float(h5.attrs["fmcw_chirp_duration"]), CHIRP_DURATION)
        assert np.isclose(float(h5.attrs["fmcw_chirp_period"]), CHIRP_PERIOD)
        assert metadata["fmcw_dechirp_mode"] == "none"
        assert metadata["fmcw_sources"][0]["transmitter_name"] == "BistaticFmcwTx"

    if fullscale <= 0.0:
        raise ValueError("simulation produced an all-zero FMCW result")
    return iq, fs, start_time, fullscale, metadata


def target_position(t):
    return TARGET_START + TARGET_VELOCITY * np.asarray(t)[:, None]


def bistatic_delay(t):
    pos = target_position(t)
    tx_leg = np.linalg.norm(pos - TX_POS, axis=1)
    rx_leg = np.linalg.norm(RX_POS - pos, axis=1)
    return (tx_leg + rx_leg) / C


def bistatic_path_rate(pos):
    u_tx = (pos - TX_POS) / np.linalg.norm(pos - TX_POS)
    u_rx = (pos - RX_POS) / np.linalg.norm(pos - RX_POS)
    return float(np.dot(TARGET_VELOCITY, u_tx + u_rx))


def slow_time_path_rate_margin():
    wavelength = C / FC
    unambiguous_path_rate = wavelength / (2.0 * CHIRP_PERIOD)
    target_end = TARGET_START + TARGET_VELOCITY * CHIRP_COUNT * CHIRP_PERIOD
    path_rate = max(abs(bistatic_path_rate(TARGET_START)), abs(bistatic_path_rate(target_end)))
    return path_rate, unambiguous_path_rate


def reference_phase(t):
    chirp_index = np.floor(t / CHIRP_PERIOD)
    u = t - chirp_index * CHIRP_PERIOD
    return 2.0 * np.pi * START_FREQUENCY_OFFSET * u + np.pi * CHIRP_RATE * u * u


def received_phase(t, tau):
    t_ret = t - tau
    chirp_index = np.floor(t_ret / CHIRP_PERIOD)
    u_ret = t_ret - chirp_index * CHIRP_PERIOD
    return (
        -2.0 * np.pi * FC * tau
        + 2.0 * np.pi * START_FREQUENCY_OFFSET * u_ret
        + np.pi * CHIRP_RATE * u_ret * u_ret
    )


def expected_if_phase(t):
    tau = bistatic_delay(t)
    return reference_phase(t) - received_phase(t, tau)


def external_dechirp(raw_iq, fs, start_time):
    times = start_time + np.arange(len(raw_iq)) / fs
    reference = np.exp(1j * reference_phase(times))
    return reference * np.conj(raw_iq)


def fit_frequency(t, values):
    phase = np.unwrap(np.angle(values)) if np.iscomplexobj(values) else np.asarray(values)
    x = t - np.mean(t)
    y = phase - np.mean(phase)
    slope = np.dot(x, y) / np.dot(x, x)
    return slope / (2.0 * np.pi)


def estimate_chirp_frequencies(if_iq, fs, start_time):
    target_end = TARGET_START + TARGET_VELOCITY * CHIRP_COUNT * CHIRP_PERIOD
    start_path = np.linalg.norm(TARGET_START - TX_POS) + np.linalg.norm(RX_POS - TARGET_START)
    end_path = np.linalg.norm(target_end - TX_POS) + np.linalg.norm(RX_POS - target_end)
    max_path = max(start_path, end_path)
    guard = max(10.0e-6, 4.0 * max_path / C)
    measured = []
    expected = []
    chirp_times = []

    for chirp in range(CHIRP_COUNT):
        t0 = chirp * CHIRP_PERIOD + guard
        t1 = chirp * CHIRP_PERIOD + CHIRP_DURATION - guard
        i0 = int(np.ceil((t0 - start_time) * fs))
        i1 = int(np.floor((t1 - start_time) * fs))
        if i1 <= i0 + 32:
            continue
        times = start_time + np.arange(i0, i1) / fs
        measured.append(fit_frequency(times, if_iq[i0:i1]))
        expected.append(fit_frequency(times, expected_if_phase(times)))
        chirp_times.append(0.5 * (t0 + t1))

    return np.array(chirp_times), np.array(measured), np.array(expected)


def estimate_slow_time_frequency(if_iq, fs, start_time):
    local_time = 0.5 * CHIRP_DURATION
    chirp_times = np.arange(CHIRP_COUNT) * CHIRP_PERIOD + local_time
    indices = np.rint((chirp_times - start_time) * fs).astype(int)
    valid = (indices >= 0) & (indices < len(if_iq))
    sample_times = start_time + indices[valid] / fs
    measured = fit_frequency(sample_times, if_iq[indices[valid]])
    expected = fit_frequency(sample_times, expected_if_phase(sample_times))
    nyquist = 0.5 / CHIRP_PERIOD
    return measured, expected, nyquist


def make_plot(output_path, chirp_times, measured, expected, raw_iq, if_iq, fs, start_time):
    fig, axes = plt.subplots(3, 1, figsize=(11, 11), constrained_layout=True)

    axes[0].plot(chirp_times * 1e3, expected, "k--", label="analytic external IF")
    axes[0].plot(chirp_times * 1e3, measured, "o", ms=4, label="measured external IF")
    axes[0].set_title("Bistatic FMCW External Dechirp")
    axes[0].set_xlabel("Time (ms)")
    axes[0].set_ylabel("Beat frequency (Hz)")
    axes[0].grid(True, alpha=0.35)
    axes[0].legend()

    chirp = CHIRP_COUNT // 2
    i0 = int((chirp * CHIRP_PERIOD - start_time) * fs)
    i1 = int(((chirp + 1) * CHIRP_PERIOD - start_time) * fs)
    window = np.hanning(i1 - i0)
    freqs = np.fft.fftshift(np.fft.fftfreq(i1 - i0, d=1.0 / fs))

    raw_spectrum = np.fft.fftshift(np.fft.fft(raw_iq[i0:i1] * window))
    raw_db = 20.0 * np.log10(np.abs(raw_spectrum) + 1.0e-30)
    raw_db -= np.max(raw_db)
    axes[1].plot(freqs / 1e3, raw_db)
    axes[1].set_xlim(-650.0, 650.0)
    axes[1].set_ylim(-90.0, 3.0)
    axes[1].set_title(f"Raw Undechirped Baseband Spectrum (chirp {chirp})")
    axes[1].set_xlabel("Frequency (kHz)")
    axes[1].set_ylabel("Relative magnitude (dB)")
    axes[1].grid(True, alpha=0.35)

    if_spectrum = np.fft.fftshift(np.fft.fft(if_iq[i0:i1] * window))
    if_db = 20.0 * np.log10(np.abs(if_spectrum) + 1.0e-30)
    if_db -= np.max(if_db)
    axes[2].plot(freqs, if_db)
    axes[2].axvline(expected[chirp], color="k", linestyle="--", label="analytic IF")
    axes[2].set_xlim(0.0, 3_000.0)
    axes[2].set_ylim(-90.0, 3.0)
    axes[2].set_title("Externally Dechirped IF Spectrum")
    axes[2].set_xlabel("Frequency (Hz)")
    axes[2].set_ylabel("Relative magnitude (dB)")
    axes[2].grid(True, alpha=0.35)
    axes[2].legend()

    fig.savefig(output_path, dpi=180)


def main():
    parser = argparse.ArgumentParser(description="Analyze the bistatic raw FMCW example with external dechirping.")
    parser.add_argument("--results-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    output_dir = (args.output_dir or results_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_iq, fs, start_time, fullscale, metadata = load_raw_iq(results_dir / RESULT_FILE)
    path_rate, unambiguous_path_rate = slow_time_path_rate_margin()
    if path_rate >= unambiguous_path_rate:
        raise SystemExit(
            f"scenario violates bistatic slow-time Nyquist: {path_rate:.3f} >= {unambiguous_path_rate:.3f} m/s"
        )

    if_iq = external_dechirp(raw_iq, fs, start_time)
    chirp_times, measured, expected = estimate_chirp_frequencies(if_iq, fs, start_time)
    slow_measured, slow_expected, slow_nyquist = estimate_slow_time_frequency(if_iq, fs, start_time)
    errors = measured - expected

    mae = float(np.mean(np.abs(errors)))
    max_error = float(np.max(np.abs(errors)))
    first_error = float(errors[0])
    last_error = float(errors[-1])

    print("Bistatic external dechirp verification")
    print(f"Samples: {len(raw_iq)} at {fs:.1f} Hz, fullscale={fullscale:.6e}")
    print(f"Metadata mode: {metadata['mode']}, dechirp={metadata['fmcw_dechirp_mode']}")
    print(f"Path-rate / unambiguous path-rate: {path_rate:.3f} / {unambiguous_path_rate:.3f} m/s")
    print(
        f"Slow-time frequency: measured={slow_measured:.2f} Hz, "
        f"expected={slow_expected:.2f} Hz, Nyquist=+/-{slow_nyquist:.2f} Hz"
    )
    print(f"Measured IF start/end: {measured[0]:.2f} Hz -> {measured[-1]:.2f} Hz")
    print(f"Expected IF start/end: {expected[0]:.2f} Hz -> {expected[-1]:.2f} Hz")
    print(f"IF error MAE/max: {mae:.3f} Hz / {max_error:.3f} Hz")
    print(f"Endpoint errors: {first_error:.3f} Hz, {last_error:.3f} Hz")

    if mae > 20.0 or max_error > 50.0:
        raise SystemExit("verification failed: externally dechirped IF does not match analytic bistatic model")
    if abs(slow_expected) >= slow_nyquist or abs(slow_measured - slow_expected) > 5.0:
        raise SystemExit("verification failed: inter-chirp phase does not match unambiguous slow-time model")

    plot_path = output_dir / PLOT_FILE
    make_plot(plot_path, chirp_times, measured, expected, raw_iq, if_iq, fs, start_time)
    print(f"Saved {plot_path}")


if __name__ == "__main__":
    main()
