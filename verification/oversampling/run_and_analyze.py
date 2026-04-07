#!/usr/bin/env python3
"""
Run the standalone oversampling audit harness across a parameter matrix and
analyze the generated CSV/manifest outputs.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
DEFAULT_BUILD_DIR = ROOT / "build"
DEFAULT_EXECUTABLE = DEFAULT_BUILD_DIR / "src" / "oversampling_audit"
DEFAULT_OUTPUT_DIR = ROOT / "output" / "analysis_runs"
DEFAULT_ANALYSIS_DIRNAME = "analysis"
DEFAULT_RATIOS = [1, 2, 4, 8, 16, 24, 32, 64]
DEFAULT_FILTER_LENGTHS = [8, 16, 33, 64]
DEFAULT_SEED = 42


@dataclass(frozen=True)
class CaseDescriptor:
    ratio: int
    filter_length: int
    suite: str
    path: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the oversampling audit executable over a matrix and analyze the outputs."
    )
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=DEFAULT_BUILD_DIR,
        help="CMake build directory for the standalone harness.",
    )
    parser.add_argument(
        "--executable",
        type=Path,
        default=DEFAULT_EXECUTABLE,
        help="Path to the built oversampling_audit executable.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where raw harness outputs will be written.",
    )
    parser.add_argument(
        "--analysis-dirname",
        default=DEFAULT_ANALYSIS_DIRNAME,
        help="Subdirectory name within output-dir for aggregated analysis artifacts.",
    )
    parser.add_argument(
        "--ratios",
        default=",".join(str(v) for v in DEFAULT_RATIOS),
        help="Comma-separated oversampling ratios to run.",
    )
    parser.add_argument(
        "--filter-lengths",
        default=",".join(str(v) for v in DEFAULT_FILTER_LENGTHS),
        help="Comma-separated render filter lengths to run.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help="Seed passed to the C++ harness.",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Do not run cmake --build before executing the harness.",
    )
    parser.add_argument(
        "--skip-run",
        action="store_true",
        help="Analyze an existing output tree without invoking the harness.",
    )
    parser.add_argument(
        "--clean-output",
        action="store_true",
        help="Delete the output directory before running.",
    )
    return parser.parse_args()


def parse_int_list(value: str) -> list[int]:
    out = [int(token.strip()) for token in value.split(",") if token.strip()]
    if not out:
        raise ValueError("expected a non-empty comma-separated integer list")
    return out


def run_command(command: list[str], cwd: Path) -> None:
    print("+", " ".join(str(part) for part in command))
    subprocess.run(command, cwd=cwd, check=True)


def build_harness(build_dir: Path) -> None:
    if not build_dir.exists():
        raise FileNotFoundError(f"Build directory does not exist: {build_dir}")
    run_command(["cmake", "--build", str(build_dir)], cwd=ROOT)


def run_matrix(executable: Path, output_dir: Path, ratios: list[int], filter_lengths: list[int], seed: int) -> None:
    if not executable.exists():
        raise FileNotFoundError(f"Executable does not exist: {executable}")
    command = [
        str(executable),
        "matrix",
        "--ratios",
        ",".join(str(v) for v in ratios),
        "--filter-lengths",
        ",".join(str(v) for v in filter_lengths),
        "--seed",
        str(seed),
        "--output-dir",
        str(output_dir),
    ]
    run_command(command, cwd=ROOT)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_summary(path: Path) -> dict[str, str]:
    rows = read_csv(path)
    return {row["metric"]: row["value"] for row in rows}


def read_manifest(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def as_float(value: str) -> float:
    return float(value)


def safe_rel_error(observed: float, expected: float) -> float:
    if expected == 0.0:
        return abs(observed - expected)
    return abs(observed - expected) / abs(expected)


def find_cases(output_dir: Path) -> list[CaseDescriptor]:
    cases: list[CaseDescriptor] = []
    for manifest_path in output_dir.glob("ratio_*/filter_*/*/manifest.json"):
        suite_dir = manifest_path.parent
        suite = suite_dir.name
        filter_token = suite_dir.parent.name
        ratio_token = suite_dir.parent.parent.name
        ratio = int(ratio_token.removeprefix("ratio_"))
        filter_length = int(filter_token.removeprefix("filter_"))
        cases.append(
            CaseDescriptor(
                ratio=ratio,
                filter_length=filter_length,
                suite=suite,
                path=suite_dir,
            )
        )
    cases.sort(key=lambda item: (item.ratio, item.filter_length, item.suite))
    return cases


def analyze_resampler(case: CaseDescriptor) -> dict[str, Any]:
    summary = read_summary(case.path / "summary.csv")
    tones = read_csv(case.path / "tone_roundtrip.csv")
    scenarios = read_csv(case.path / "scenario_roundtrip.csv")
    alignment = read_csv(case.path / "group_delay_alignment.csv")

    min_tone_corr = min(as_float(row["normalized_correlation"]) for row in tones)
    max_tone_rms = max(as_float(row["rms_error"]) for row in tones)
    min_scenario_corr = min(as_float(row["normalized_correlation"]) for row in scenarios)
    max_scenario_rms = max(as_float(row["rms_error"]) for row in scenarios)
    max_peak_error = max(abs(as_float(row["peak_index_error"])) for row in alignment)
    fir_sum = as_float(summary["fir_sum"])

    return {
        "suite": case.suite,
        "ratio": case.ratio,
        "filter_length": case.filter_length,
        "fir_sum_error": fir_sum - case.ratio,
        "min_tone_correlation": min_tone_corr,
        "max_tone_rms_error": max_tone_rms,
        "min_scenario_correlation": min_scenario_corr,
        "max_scenario_rms_error": max_scenario_rms,
        "max_group_delay_peak_error": max_peak_error,
        "pass": (
            abs(fir_sum - case.ratio) < 2e-3
            and min_tone_corr > 0.98
            and min_scenario_corr > 0.95
            and max_peak_error <= 0.0
        ),
    }


def analyze_render(case: CaseDescriptor) -> dict[str, Any]:
    summary = read_summary(case.path / "summary.csv")
    constant_rows = read_csv(case.path / "constant_render.csv")
    interpolation_rows = read_csv(case.path / "power_phase_interpolation.csv")
    fractional_rows = read_csv(case.path / "fractional_delay_sweep.csv")
    oversampled_rows = read_csv(case.path / "oversampled_render_comparison.csv")

    center_index = case.filter_length
    center_row = next(row for row in constant_rows if int(row["index"]) == center_index)
    constant_real_error = abs(as_float(center_row["real"]) - as_float(summary["constant_expected_real"]))
    constant_imag_error = abs(as_float(center_row["imag"]) - as_float(summary["constant_expected_imag"]))

    max_interp_error = 0.0
    for row in interpolation_rows:
        real_error = abs(as_float(row["observed_real"]) - as_float(row["expected_real"]))
        imag_error = abs(as_float(row["observed_imag"]) - as_float(row["expected_imag"]))
        max_interp_error = max(max_interp_error, real_error, imag_error)

    max_fractional_error = 0.0
    for row in fractional_rows:
        max_fractional_error = max(
            max_fractional_error,
            abs(as_float(row["observed_real"]) - as_float(row["expected_real"])),
            abs(as_float(row["observed_imag"])),
        )

    oversampled_corr = min(as_float(row["normalized_correlation"]) for row in oversampled_rows)
    oversampled_rms = max(as_float(row["rms_error"]) for row in oversampled_rows)

    return {
        "suite": case.suite,
        "ratio": case.ratio,
        "filter_length": case.filter_length,
        "constant_center_real_error": constant_real_error,
        "constant_center_imag_error": constant_imag_error,
        "max_interpolation_error": max_interp_error,
        "max_fractional_delay_error": max_fractional_error,
        "oversampled_render_correlation": oversampled_corr,
        "oversampled_render_rms_error": oversampled_rms,
        "pass": (
            constant_real_error < 1e-5
            and constant_imag_error < 1e-5
            and max_fractional_error < 1e-3
            and oversampled_corr > 0.999
            and oversampled_rms < 5e-2
        ),
    }


def analyze_pipeline(case: CaseDescriptor) -> dict[str, Any]:
    summary = read_summary(case.path / "summary.csv")
    thermal_rows = read_csv(case.path / "thermal_noise_stats.csv")
    downsampling_rows = read_csv(case.path / "downsampling_quantization.csv")

    max_thermal_rel_error = 0.0
    for row in thermal_rows:
        expected = as_float(row["expected_per_channel_power"])
        variance_real = as_float(row["variance_real"])
        variance_imag = as_float(row["variance_imag"])
        if expected == 0.0:
            continue
        max_thermal_rel_error = max(
            max_thermal_rel_error,
            safe_rel_error(variance_real, expected),
            safe_rel_error(variance_imag, expected),
        )

    min_downsampling_corr = min(as_float(row["normalized_correlation"]) for row in downsampling_rows)
    max_downsampling_rms = max(as_float(row["rms_error"]) for row in downsampling_rows)

    histogram_level_count: int | None = None
    histogram_paths = sorted(case.path.glob("adc_histogram_bits_*.csv"))
    if histogram_paths:
        histogram_level_count = max(len(read_csv(path)) for path in histogram_paths)

    return {
        "suite": case.suite,
        "ratio": case.ratio,
        "filter_length": case.filter_length,
        "thermal_max_relative_error": max_thermal_rel_error,
        "downsampling_correlation": min_downsampling_corr,
        "downsampling_rms_error": max_downsampling_rms,
        "downsampled_fullscale": as_float(summary["downsampled_fullscale"]),
        "max_adc_histogram_levels": histogram_level_count if histogram_level_count is not None else 0,
        "pass": max_thermal_rel_error < 0.25 and min_downsampling_corr > 0.99 and max_downsampling_rms < 0.04,
    }


def analyze_clocking(case: CaseDescriptor) -> dict[str, Any]:
    summary = read_summary(case.path / "summary.csv")
    tx_rows = read_csv(case.path / "transmitter_prf_quantization.csv")
    rx_rows = read_csv(case.path / "receiver_window_quantization.csv")
    cw_rows = read_csv(case.path / "cw_timestep_and_buffer.csv")

    worst_tx = max(tx_rows, key=lambda row: abs(as_float(row["error"])))
    worst_rx = max(rx_rows, key=lambda row: abs(as_float(row["skip_error"])))
    smallest_dt = min(as_float(row["dt_sim"]) for row in cw_rows)
    largest_buffer = max(int(row["cw_buffer_size"]) for row in cw_rows)

    return {
        "suite": case.suite,
        "ratio": case.ratio,
        "filter_length": case.filter_length,
        "max_abs_prf_error": abs(as_float(summary["max_abs_prf_error"])),
        "max_abs_skip_error": abs(as_float(summary["max_abs_skip_error"])),
        "worst_prf_requested": as_float(worst_tx["requested_prf"]),
        "worst_skip_requested": as_float(worst_rx["requested_skip"]),
        "smallest_dt_sim": smallest_dt,
        "largest_cw_buffer_size": largest_buffer,
        "pass": True,
    }


SUITE_ANALYZERS = {
    "resampler": analyze_resampler,
    "render": analyze_render,
    "pipeline": analyze_pipeline,
    "clocking": analyze_clocking,
}


def analyze_cases(cases: list[CaseDescriptor]) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for case in cases:
        manifest = read_manifest(case.path / "manifest.json")
        if manifest["suite"] != case.suite:
            raise RuntimeError(f"Manifest suite mismatch at {case.path}")
        analyzer = SUITE_ANALYZERS.get(case.suite)
        if analyzer is None:
            raise RuntimeError(f"No analyzer registered for suite {case.suite}")
        results.append(analyzer(case))
    return results


def coerce_table_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.12e}"
    return str(value)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError(f"No rows to write to {path}")
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: coerce_table_value(row.get(key, "")) for key in fieldnames})


def summarize_by_suite(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    suite_names = sorted({str(row["suite"]) for row in rows})
    out: list[dict[str, Any]] = []
    for suite in suite_names:
        suite_rows = [row for row in rows if row["suite"] == suite]
        pass_count = sum(1 for row in suite_rows if bool(row["pass"]))
        out.append(
            {
                "suite": suite,
                "cases": len(suite_rows),
                "passes": pass_count,
                "failures": len(suite_rows) - pass_count,
                "pass_rate": pass_count / len(suite_rows) if suite_rows else 0.0,
            }
        )
    return out


def build_global_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    suite_summary = summarize_by_suite(rows)
    total_cases = len(rows)
    total_passes = sum(1 for row in rows if bool(row["pass"]))
    failing_cases = [row for row in rows if not bool(row["pass"])]
    return {
        "total_cases": total_cases,
        "total_passes": total_passes,
        "total_failures": total_cases - total_passes,
        "overall_pass_rate": total_passes / total_cases if total_cases else 0.0,
        "suite_summary": suite_summary,
        "failing_cases": failing_cases,
    }


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def worst_case(rows: list[dict[str, Any]], suite: str, metric: str, reverse: bool = True) -> dict[str, Any] | None:
    suite_rows = [row for row in rows if row["suite"] == suite and metric in row]
    if not suite_rows:
        return None
    return sorted(suite_rows, key=lambda row: float(row[metric]), reverse=reverse)[0]


def generate_markdown_report(rows: list[dict[str, Any]], global_summary: dict[str, Any], args: argparse.Namespace) -> str:
    lines: list[str] = []
    lines.append("# Oversampling Audit Analysis")
    lines.append("")
    lines.append("## Run Configuration")
    lines.append("")
    lines.append(f"- Output root: `{args.output_dir}`")
    lines.append(f"- Ratios: `{args.ratios}`")
    lines.append(f"- Filter lengths: `{args.filter_lengths}`")
    lines.append(f"- Seed: `{args.seed}`")
    lines.append("")
    lines.append("## Overall Result")
    lines.append("")
    lines.append(f"- Cases analyzed: `{global_summary['total_cases']}`")
    lines.append(f"- Passing cases: `{global_summary['total_passes']}`")
    lines.append(f"- Failing cases: `{global_summary['total_failures']}`")
    lines.append(f"- Overall pass rate: `{global_summary['overall_pass_rate']:.2%}`")
    lines.append("")
    lines.append("## Suite Summary")
    lines.append("")
    lines.append("| Suite | Cases | Passes | Failures | Pass Rate |")
    lines.append("| --- | ---: | ---: | ---: | ---: |")
    for item in global_summary["suite_summary"]:
        lines.append(
            f"| {item['suite']} | {item['cases']} | {item['passes']} | {item['failures']} | {item['pass_rate']:.2%} |"
        )
    lines.append("")
    lines.append("## Key Findings")
    lines.append("")

    resampler = worst_case(rows, "resampler", "max_tone_rms_error", reverse=True)
    if resampler is not None:
        lines.append(
            "- Resampler worst tone RMS error: "
            f"`{resampler['max_tone_rms_error']:.6e}` at ratio `{resampler['ratio']}`, filter `{resampler['filter_length']}`."
        )
    render = worst_case(rows, "render", "max_fractional_delay_error", reverse=True)
    if render is not None:
        lines.append(
            "- Render worst fractional-delay error: "
            f"`{render['max_fractional_delay_error']:.6e}` at ratio `{render['ratio']}`, filter `{render['filter_length']}`."
        )
    pipeline = worst_case(rows, "pipeline", "thermal_max_relative_error", reverse=True)
    if pipeline is not None:
        lines.append(
            "- Pipeline worst thermal-noise relative error: "
            f"`{pipeline['thermal_max_relative_error']:.6e}` at ratio `{pipeline['ratio']}`, filter `{pipeline['filter_length']}`."
        )
    clocking = worst_case(rows, "clocking", "max_abs_prf_error", reverse=True)
    if clocking is not None:
        lines.append(
            "- Clocking worst PRF quantization error: "
            f"`{clocking['max_abs_prf_error']:.6e}` at ratio `{clocking['ratio']}`, filter `{clocking['filter_length']}`."
        )
    lines.append("")

    if global_summary["failing_cases"]:
        lines.append("## Failing Cases")
        lines.append("")
        lines.append("| Suite | Ratio | Filter | Notes |")
        lines.append("| --- | ---: | ---: | --- |")
        for row in global_summary["failing_cases"]:
            notes = []
            if row["suite"] == "resampler":
                notes.append(f"min tone corr={row['min_tone_correlation']:.6e}")
                notes.append(f"min scenario corr={row['min_scenario_correlation']:.6e}")
            elif row["suite"] == "render":
                notes.append(f"max fractional error={row['max_fractional_delay_error']:.6e}")
                notes.append(f"oversampled corr={row['oversampled_render_correlation']:.6e}")
            elif row["suite"] == "pipeline":
                notes.append(f"thermal rel err={row['thermal_max_relative_error']:.6e}")
                notes.append(f"downsample corr={row['downsampling_correlation']:.6e}")
            lines.append(f"| {row['suite']} | {row['ratio']} | {row['filter_length']} | {'; '.join(notes)} |")
        lines.append("")
    else:
        lines.append("## Failing Cases")
        lines.append("")
        lines.append("No failing cases were detected using the current analysis thresholds.")
        lines.append("")

    lines.append("## Analysis Artifacts")
    lines.append("")
    lines.append(f"- `analysis/case_analysis.csv`: per-case extracted metrics and pass/fail results")
    lines.append(f"- `analysis/suite_summary.csv`: aggregated pass/fail counts by suite")
    lines.append(f"- `analysis/global_summary.json`: machine-readable rollup of all analysis results")
    lines.append(f"- `analysis/report.md`: this report")
    lines.append("")
    return "\n".join(lines)


def write_analysis_outputs(rows: list[dict[str, Any]], args: argparse.Namespace) -> None:
    analysis_dir = args.output_dir / args.analysis_dirname
    analysis_dir.mkdir(parents=True, exist_ok=True)

    case_analysis_path = analysis_dir / "case_analysis.csv"
    suite_summary_path = analysis_dir / "suite_summary.csv"
    global_summary_path = analysis_dir / "global_summary.json"
    report_path = analysis_dir / "report.md"

    suite_summary = summarize_by_suite(rows)
    global_summary = build_global_summary(rows)

    write_csv(case_analysis_path, rows)
    write_csv(suite_summary_path, suite_summary)
    write_json(global_summary_path, global_summary)
    report_text = generate_markdown_report(rows, global_summary, args)
    report_path.write_text(report_text + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.build_dir = args.build_dir.resolve()
    args.executable = args.executable.resolve()
    args.output_dir = args.output_dir.resolve()
    ratios = parse_int_list(args.ratios)
    filter_lengths = parse_int_list(args.filter_lengths)

    if args.clean_output and args.output_dir.exists():
        shutil.rmtree(args.output_dir)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_build:
        build_harness(args.build_dir)

    if not args.skip_run:
        run_matrix(args.executable, args.output_dir, ratios, filter_lengths, args.seed)

    cases = find_cases(args.output_dir)
    if not cases:
        raise RuntimeError(f"No matrix cases found under {args.output_dir}")

    analysis_rows = analyze_cases(cases)
    write_analysis_outputs(analysis_rows, args)

    global_summary = build_global_summary(analysis_rows)
    print(
        "Analysis complete:",
        f"{global_summary['total_passes']}/{global_summary['total_cases']} cases passed",
        f"({global_summary['overall_pass_rate']:.2%})",
    )
    print(f"Artifacts written under: {args.output_dir / args.analysis_dirname}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
