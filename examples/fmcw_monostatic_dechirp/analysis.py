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
RESULT_FILE = "MonoFmcwRadar_results.h5"
PLOT_FILE = "monostatic_dechirp_analysis.png"

FC = 10.0e9
CHIRP_BANDWIDTH = 1.0e6
CHIRP_DURATION = 1.0e-3
CHIRP_PERIOD = 1.0e-3
CHIRP_COUNT = 64
START_FREQUENCY_OFFSET = -5.0e5
CHIRP_RATE = CHIRP_BANDWIDTH / CHIRP_DURATION

RADAR_POS = np.array([0.0, 0.0, 0.0])
TARGET_START = np.array([150.0, 0.0, 0.0])
TARGET_VELOCITY = np.array([5.0, 0.0, 0.0])


def attr_text(attrs, name):
    value = attrs[name]
    return value.decode("utf-8") if isinstance(value, bytes) else value


def load_iq(path):
    if not path.exists():
        raise FileNotFoundError(f"missing simulation output: {path}")

    with h5py.File(path, "r") as h5:
        fs = float(h5.attrs["sampling_rate"])
        start_time = float(h5.attrs["start_time"])
        fullscale = float(h5.attrs["fullscale"])
        metadata = json.loads(attr_text(h5.attrs, "fers_metadata_json"))
        iq = (h5["I_data"][:] + 1j * h5["Q_data"][:]) * fullscale

        assert attr_text(h5.attrs, "data_mode") == "fmcw"
        assert attr_text(h5.attrs, "fmcw_dechirp_mode") == "physical"
        assert attr_text(h5.attrs, "fmcw_dechirp_reference_source") == "attached"
        assert attr_text(h5.attrs, "fmcw_dechirp_reference_transmitter_name") == "MonoFmcwRadar"
        assert int(h5.attrs["fmcw_source_count"]) == 1
        assert int(h5.attrs["fmcw_chirp_count"]) == CHIRP_COUNT
        assert np.isclose(float(h5.attrs["fmcw_chirp_bandwidth"]), CHIRP_BANDWIDTH)
        assert np.isclose(float(h5.attrs["fmcw_chirp_duration"]), CHIRP_DURATION)
        assert np.isclose(float(h5.attrs["fmcw_chirp_period"]), CHIRP_PERIOD)
        assert metadata["fmcw_dechirp_mode"] == "physical"
        assert metadata["fmcw_dechirp_reference_source"] == "attached"

    if fullscale <= 0.0:
        raise ValueError("simulation produced an all-zero FMCW result")
    return iq, fs, start_time, fullscale, metadata


def target_position(t):
    return TARGET_START + TARGET_VELOCITY * np.asarray(t)[:, None]


def monostatic_delay(t):
    ranges = np.linalg.norm(target_position(t) - RADAR_POS, axis=1)
    return 2.0 * ranges / C


def slow_time_velocity_margin():
    wavelength = C / FC
    unambiguous_velocity = wavelength / (4.0 * CHIRP_PERIOD)
    radial_unit = (TARGET_START - RADAR_POS) / np.linalg.norm(TARGET_START - RADAR_POS)
    radial_speed = abs(float(np.dot(TARGET_VELOCITY, radial_unit)))
    return radial_speed, unambiguous_velocity


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
    tau = monostatic_delay(t)
    return reference_phase(t) - received_phase(t, tau)


def fit_frequency(t, values):
    phase = np.unwrap(np.angle(values)) if np.iscomplexobj(values) else np.asarray(values)
    x = t - np.mean(t)
    y = phase - np.mean(phase)
    slope = np.dot(x, y) / np.dot(x, x)
    return slope / (2.0 * np.pi)


def estimate_chirp_frequencies(iq, fs, start_time):
    max_tau = 2.0 * np.linalg.norm(TARGET_START + TARGET_VELOCITY * CHIRP_COUNT * CHIRP_PERIOD) / C
    guard = max(8.0e-6, 4.0 * max_tau)
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
        measured.append(fit_frequency(times, iq[i0:i1]))
        expected.append(fit_frequency(times, expected_if_phase(times)))
        chirp_times.append(0.5 * (t0 + t1))

    return np.array(chirp_times), np.array(measured), np.array(expected)


def estimate_slow_time_frequency(iq, fs, start_time):
    local_time = 0.5 * CHIRP_DURATION
    chirp_times = np.arange(CHIRP_COUNT) * CHIRP_PERIOD + local_time
    indices = np.rint((chirp_times - start_time) * fs).astype(int)
    valid = (indices >= 0) & (indices < len(iq))
    sample_times = start_time + indices[valid] / fs
    measured = fit_frequency(sample_times, iq[indices[valid]])
    expected = fit_frequency(sample_times, expected_if_phase(sample_times))
    nyquist = 0.5 / CHIRP_PERIOD
    return measured, expected, nyquist


def make_plot(output_path, chirp_times, measured, expected, iq, fs, start_time):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), constrained_layout=True)

    axes[0].plot(chirp_times * 1e3, expected, "k--", label="analytic IF")
    axes[0].plot(chirp_times * 1e3, measured, "o", ms=4, label="measured IF")
    axes[0].set_title("Monostatic FMCW Physical Dechirp")
    axes[0].set_xlabel("Time (ms)")
    axes[0].set_ylabel("Beat frequency (Hz)")
    axes[0].grid(True, alpha=0.35)
    axes[0].legend()

    chirp = CHIRP_COUNT // 2
    i0 = int((chirp * CHIRP_PERIOD - start_time) * fs)
    i1 = int(((chirp + 1) * CHIRP_PERIOD - start_time) * fs)
    segment = iq[i0:i1]
    freqs = np.fft.fftshift(np.fft.fftfreq(len(segment), d=1.0 / fs))
    spectrum = np.fft.fftshift(np.fft.fft(segment * np.hanning(len(segment)), n=len(segment)))
    spectrum_db = 20.0 * np.log10(np.abs(spectrum) + 1.0e-30)
    spectrum_db -= np.max(spectrum_db)
    axes[1].plot(freqs, spectrum_db)
    axes[1].axvline(expected[chirp], color="k", linestyle="--", label="analytic IF")
    axes[1].set_xlim(0.0, 5_000.0)
    axes[1].set_ylim(-90.0, 3.0)
    axes[1].set_title(f"Mid-Chirp IF Spectrum (chirp {chirp})")
    axes[1].set_xlabel("Frequency (Hz)")
    axes[1].set_ylabel("Relative magnitude (dB)")
    axes[1].grid(True, alpha=0.35)
    axes[1].legend()

    fig.savefig(output_path, dpi=180)


def main():
    parser = argparse.ArgumentParser(description="Analyze the monostatic built-in FMCW dechirp example.")
    parser.add_argument("--results-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    output_dir = (args.output_dir or results_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    iq, fs, start_time, fullscale, metadata = load_iq(results_dir / RESULT_FILE)
    radial_speed, unambiguous_velocity = slow_time_velocity_margin()
    if radial_speed >= unambiguous_velocity:
        raise SystemExit(
            f"scenario violates monostatic slow-time Nyquist: {radial_speed:.3f} >= {unambiguous_velocity:.3f} m/s"
        )

    chirp_times, measured, expected = estimate_chirp_frequencies(iq, fs, start_time)
    slow_measured, slow_expected, slow_nyquist = estimate_slow_time_frequency(iq, fs, start_time)
    errors = measured - expected

    mae = float(np.mean(np.abs(errors)))
    max_error = float(np.max(np.abs(errors)))
    first_error = float(errors[0])
    last_error = float(errors[-1])

    print("Monostatic physical dechirp verification")
    print(f"Samples: {len(iq)} at {fs:.1f} Hz, fullscale={fullscale:.6e}")
    print(f"Metadata mode: {metadata['mode']}, dechirp={metadata['fmcw_dechirp_mode']}")
    print(f"Radial speed / unambiguous speed: {radial_speed:.3f} / {unambiguous_velocity:.3f} m/s")
    print(
        f"Slow-time frequency: measured={slow_measured:.2f} Hz, "
        f"expected={slow_expected:.2f} Hz, Nyquist=+/-{slow_nyquist:.2f} Hz"
    )
    print(f"Measured IF start/end: {measured[0]:.2f} Hz -> {measured[-1]:.2f} Hz")
    print(f"Expected IF start/end: {expected[0]:.2f} Hz -> {expected[-1]:.2f} Hz")
    print(f"IF error MAE/max: {mae:.3f} Hz / {max_error:.3f} Hz")
    print(f"Endpoint errors: {first_error:.3f} Hz, {last_error:.3f} Hz")

    if mae > 15.0 or max_error > 40.0:
        raise SystemExit("verification failed: measured IF does not match analytic physical dechirp model")
    if abs(slow_expected) >= slow_nyquist or abs(slow_measured - slow_expected) > 5.0:
        raise SystemExit("verification failed: inter-chirp phase does not match unambiguous slow-time model")

    plot_path = output_dir / PLOT_FILE
    make_plot(plot_path, chirp_times, measured, expected, iq, fs, start_time)
    print(f"Saved {plot_path}")


if __name__ == "__main__":
    main()
