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
RESULT_FILE = "SfcwRadar_results.h5"
PLOT_FILE = "sfcw_range_profile.png"
SUMMARY_FILE = "sfcw_analysis_summary.json"

TARGET_RANGE_M = 100.0
STEP_COUNT = 16
STEP_SIZE_HZ = 100_000.0
DWELL_TIME_S = 1.0e-4
STEP_PERIOD_S = 2.0e-4
SWEEP_COUNT = 4


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

        assert attr_text(h5.attrs, "data_mode") == "sfcw"
        assert int(h5.attrs["sfcw_step_count"]) == STEP_COUNT
        assert int(h5.attrs["sfcw_sweep_count"]) == SWEEP_COUNT
        assert int(h5.attrs["sfcw_source_count"]) == 1
        assert np.isclose(float(h5.attrs["sfcw_step_size"]), STEP_SIZE_HZ)
        assert np.isclose(float(h5.attrs["sfcw_dwell_time"]), DWELL_TIME_S)
        assert np.isclose(float(h5.attrs["sfcw_step_period"]), STEP_PERIOD_S)
        assert metadata["mode"] == "sfcw"
        assert metadata["sfcw"]["step_count"] == STEP_COUNT
        assert metadata["sfcw"]["sweep_count"] == SWEEP_COUNT

    if fullscale <= 0.0:
        raise ValueError("simulation produced an all-zero SFCW result")
    return iq, fs, start_time, fullscale, metadata


def extract_step_samples(iq, fs, start_time):
    delay = 2.0 * TARGET_RANGE_M / C
    samples = np.zeros((SWEEP_COUNT, STEP_COUNT), dtype=np.complex128)
    sweep_period = STEP_COUNT * STEP_PERIOD_S

    for sweep in range(SWEEP_COUNT):
        for step in range(STEP_COUNT):
            t = sweep * sweep_period + step * STEP_PERIOD_S + delay + 0.5 * DWELL_TIME_S
            index = int(round((t - start_time) * fs))
            if index < 0 or index >= len(iq):
                raise IndexError(f"step sample outside result: sweep={sweep}, step={step}, index={index}")
            samples[sweep, step] = iq[index]

    return samples


def range_profile(step_response, nfft=1024):
    profile = np.fft.ifft(step_response, n=nfft)
    ranges = np.arange(nfft) * C / (2.0 * nfft * abs(STEP_SIZE_HZ))
    magnitude_db = 20.0 * np.log10(np.abs(profile) + 1.0e-30)
    magnitude_db -= float(np.max(magnitude_db))
    return ranges, magnitude_db


def make_plot(output_path, ranges, magnitude_db, estimated_range):
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    ax.plot(ranges, magnitude_db)
    ax.axvline(TARGET_RANGE_M, color="k", linestyle="--", label="target range")
    ax.axvline(estimated_range, color="tab:red", linestyle=":", label="SFCW peak")
    ax.set_xlim(0.0, 400.0)
    ax.set_ylim(-60.0, 3.0)
    ax.set_title("SFCW Range Profile")
    ax.set_xlabel("Range (m)")
    ax.set_ylabel("Relative magnitude (dB)")
    ax.grid(True, alpha=0.35)
    ax.legend()
    fig.savefig(output_path, dpi=180)


def main():
    parser = argparse.ArgumentParser(description="Analyze the monostatic SFCW example.")
    parser.add_argument("--results-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    output_dir = (args.output_dir or results_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    iq, fs, start_time, fullscale, metadata = load_iq(results_dir / RESULT_FILE)
    samples = extract_step_samples(iq, fs, start_time)
    step_response = np.mean(samples, axis=0)
    ranges, magnitude_db = range_profile(step_response)

    peak_index = int(np.argmax(magnitude_db))
    estimated_range = float(ranges[peak_index])
    range_error = abs(estimated_range - TARGET_RANGE_M)
    resolution = float(metadata["sfcw"]["range_resolution"])
    unambiguous = float(metadata["sfcw"]["unambiguous_range"])

    print("Monostatic SFCW verification")
    print(f"Samples: {len(iq)} at {fs:.1f} Hz, fullscale={fullscale:.6e}")
    print(f"Metadata mode: {metadata['mode']}")
    print(f"Resolution / unambiguous range: {resolution:.3f} / {unambiguous:.3f} m")
    print(f"Estimated target range: {estimated_range:.3f} m")
    print(f"Range error: {range_error:.3f} m")

    if range_error > 0.5 * resolution:
        raise SystemExit("verification failed: SFCW range profile peak is not near the target")

    plot_path = output_dir / PLOT_FILE
    make_plot(plot_path, ranges, magnitude_db, estimated_range)
    print(f"Saved {plot_path}")

    summary = {
        "mode": metadata["mode"],
        "samples": int(len(iq)),
        "sample_rate_hz": fs,
        "target_range_m": TARGET_RANGE_M,
        "estimated_range_m": estimated_range,
        "range_error_m": range_error,
        "range_resolution_m": resolution,
        "unambiguous_range_m": unambiguous,
        "effective_bandwidth_hz": float(metadata["sfcw"]["effective_bandwidth"]),
    }
    summary_path = output_dir / SUMMARY_FILE
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"Saved {summary_path}")


if __name__ == "__main__":
    main()
