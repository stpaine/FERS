#!/usr/bin/env python3
"""
summarize_build_log.py

Summarize a CMake/Make/Ninja C/C++ build log into a readable report.

The script can:
- parse compiler warnings, errors, and notes from a build log
- detect build targets and link steps
- count raw emitted warnings and unique warning sites
- distinguish header-origin warnings from source-file warnings
- optionally enrich the analysis with compile_commands.json
- emit plain text, Markdown, or JSON output

Examples:
    python summarize_build_log.py build.log
    python summarize_build_log.py build.log --markdown
    python summarize_build_log.py build.log --json
    python summarize_build_log.py build.log --compile-db build/release/compile_commands.json
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional


ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")

WARNING_RE = re.compile(
    r"""
    ^
    (?P<file>.*?)
    :
    (?P<line>\d+)
    :
    (?P<col>\d+)
    :
    \s+
    warning:
    \s+
    (?P<message>.*?)
    (?:\s+\[(?P<flag>[^\]]+)\])?
    \s*$
    """,
    re.VERBOSE,
)

ERROR_RE = re.compile(
    r"""
    ^
    (?P<file>.*?)
    :
    (?P<line>\d+)
    :
    (?P<col>\d+)
    :
    \s+
    error:
    \s+
    (?P<message>.*?)
    (?:\s+\[(?P<flag>[^\]]+)\])?
    \s*$
    """,
    re.VERBOSE,
)

NOTE_RE = re.compile(
    r"""
    ^
    (?P<file>.*?)
    :
    (?P<line>\d+)
    :
    (?P<col>\d+)
    :
    \s+
    note:
    \s+
    (?P<message>.*?)
    \s*$
    """,
    re.VERBOSE,
)

IN_FUNCTION_RE = re.compile(
    r"^\s*(?P<file>.*?): In (?:function|member function|constructor|destructor|lambda function)\s+[‘'].*$"
)
INSTANTIATION_RE = re.compile(r"^\s*(?P<file>.*?): In instantiation of [‘'].*$")
BUILDING_OBJ_RE = re.compile(r"^\[\s*\d+%\] Building (?:C|CXX) object (.+)$")
LINKING_RE = re.compile(r"^\[\s*\d+%\] Linking (.+)$")
BUILT_TARGET_RE = re.compile(r"^\[\s*\d+%\] Built target (.+)$")

COMPILER_CMD_RE = re.compile(r"(^|\s)(?:/usr/bin/)?(?:c\+\+|g\+\+|clang\+\+|cc|gcc|clang)(\s|$)")
WARNING_FLAG_IN_CMD_RE = re.compile(r"(?<!\S)(-W[a-zA-Z0-9][\w=-]*)(?!\S)")
STD_RE = re.compile(r"(?<!\S)-std=([^\s]+)")
OPT_RE = re.compile(r"(?<!\S)(-O[0-3sgz]|-Ofast)(?!\S)")
DEFINE_RE = re.compile(r"(?<!\S)-D([A-Za-z_]\w*)(?:=([^\s]+))?")
INCLUDE_RE = re.compile(r"(?<!\S)(?:-I|-isystem)\s*([^\s]+)")
MAKE_ERROR_RE = re.compile(r"(?:gmake|make|ninja)(?:\[\d+\])?: \*\*\*")
FATAL_RE = re.compile(r"\bfatal error\b", re.IGNORECASE)
FAILED_RE = re.compile(r"\b(?:FAILED|failed)\b")
ENTER_DIR_RE = re.compile(r"Entering directory '([^']+)'")
SOURCE_PATH_RE = re.compile(r"(?<!\S)(/[^ ]+\.(?:c|cc|cpp|cxx|c\+\+|C))(?!\S)")
OUTPUT_OBJ_RE = re.compile(r"(?<!\S)-o\s+([^\s]+(?:\.o|\.obj))(?!\S)")


@dataclass
class Diagnostic:
    kind: str
    file: str
    line: int
    col: int
    message: str
    flag: Optional[str] = None
    context: Optional[str] = None
    translation_unit: Optional[str] = None
    target_hint: Optional[str] = None

    def site_key(self) -> tuple:
        return (self.kind, self.file, self.line, self.col, self.message, self.flag)

    def broad_key(self) -> tuple:
        return (self.kind, self.file, self.line, self.message, self.flag)


@dataclass
class CompileCommand:
    directory: Optional[str]
    file: Optional[str]
    output: Optional[str]
    command: Optional[str]
    arguments: list[str] = field(default_factory=list)
    std: Optional[str] = None
    optimization: Optional[str] = None
    warning_flags: list[str] = field(default_factory=list)
    defines: list[str] = field(default_factory=list)
    includes: list[str] = field(default_factory=list)
    target_hint: Optional[str] = None


@dataclass
class BuildSummary:
    build_dir: Optional[str] = None
    build_command: Optional[str] = None
    success: Optional[bool] = None

    warnings: int = 0
    errors: int = 0
    notes: int = 0

    diagnostics: list[Diagnostic] = field(default_factory=list)
    compile_commands: list[CompileCommand] = field(default_factory=list)

    warnings_by_flag: Counter = field(default_factory=Counter)
    warnings_by_file: Counter = field(default_factory=Counter)
    errors_by_file: Counter = field(default_factory=Counter)

    built_targets: list[str] = field(default_factory=list)
    building_steps: list[str] = field(default_factory=list)
    linking_steps: list[str] = field(default_factory=list)
    directories_entered: list[str] = field(default_factory=list)

    compiler_warning_flags_seen: Counter = field(default_factory=Counter)
    standards_seen: Counter = field(default_factory=Counter)
    optimization_levels_seen: Counter = field(default_factory=Counter)

    raw_error_lines: list[str] = field(default_factory=list)
    fatal_lines: list[str] = field(default_factory=list)


def strip_ansi(text: str) -> str:
    return ANSI_RE.sub("", text)


def compact_message(msg: str, limit: int = 150) -> str:
    msg = " ".join(msg.split())
    return msg if len(msg) <= limit else msg[: limit - 1] + "…"


def shorten_path(path: str, build_dir: Optional[str]) -> str:
    p = path.strip()
    if not p:
        return p

    home = str(Path.home())
    if p.startswith(home):
        p = "~" + p[len(home):]

    if build_dir:
        bd = build_dir.rstrip("/")
        if p.startswith(bd):
            rel = p[len(bd):].lstrip("/")
            return f"<build>/{rel}" if rel else "<build>"

    return p


def normalize_path(path: Optional[str], base_dir: Optional[str] = None) -> Optional[str]:
    if not path:
        return None
    p = Path(path)
    if not p.is_absolute() and base_dir:
        p = Path(base_dir) / p
    try:
        return str(p.resolve(strict=False))
    except Exception:
        return str(p)


def extract_target_hint(output: Optional[str], directory: Optional[str]) -> Optional[str]:
    candidate = output or directory
    if not candidate:
        return None

    parts = Path(candidate).parts
    if "CMakeFiles" in parts:
        idx = parts.index("CMakeFiles")
        if idx + 1 < len(parts):
            name = parts[idx + 1]
            if name.endswith(".dir"):
                return name[:-4]

    for part in reversed(parts):
        if part.endswith(".dir"):
            return part[:-4]
    return None


def parse_compile_command_from_line(line: str) -> Optional[CompileCommand]:
    if not COMPILER_CMD_RE.search(line):
        return None
    if " -c " not in line:
        return None

    warning_flags = sorted(set(WARNING_FLAG_IN_CMD_RE.findall(line)))
    std_match = STD_RE.search(line)
    opt_match = OPT_RE.search(line)
    defines = [m.group(1) for m in DEFINE_RE.finditer(line)]
    includes = [m.group(1) for m in INCLUDE_RE.finditer(line)]
    source = None
    candidates = SOURCE_PATH_RE.findall(line)
    if candidates:
        source = candidates[-1]
    output_match = OUTPUT_OBJ_RE.search(line)
    output = output_match.group(1) if output_match else None

    return CompileCommand(
        directory=None,
        file=source,
        output=output,
        command=line,
        arguments=[],
        std=std_match.group(1) if std_match else None,
        optimization=opt_match.group(1) if opt_match else None,
        warning_flags=warning_flags,
        defines=defines,
        includes=includes,
        target_hint=extract_target_hint(output, None),
    )


def load_compile_database(path: str) -> list[CompileCommand]:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        raw = json.load(f)

    commands: list[CompileCommand] = []
    for entry in raw:
        directory = entry.get("directory")
        file_ = normalize_path(entry.get("file"), directory)
        output = normalize_path(entry.get("output"), directory) if entry.get("output") else None

        command_text = entry.get("command")
        arguments = entry.get("arguments", [])

        text_for_scan = command_text or " ".join(arguments)

        warning_flags = sorted(set(WARNING_FLAG_IN_CMD_RE.findall(text_for_scan)))
        std_match = STD_RE.search(text_for_scan)
        opt_match = OPT_RE.search(text_for_scan)
        defines = [m.group(1) for m in DEFINE_RE.finditer(text_for_scan)]
        includes = [m.group(1) for m in INCLUDE_RE.finditer(text_for_scan)]

        commands.append(
            CompileCommand(
                directory=normalize_path(directory),
                file=file_,
                output=output,
                command=command_text,
                arguments=arguments,
                std=std_match.group(1) if std_match else None,
                optimization=opt_match.group(1) if opt_match else None,
                warning_flags=warning_flags,
                defines=defines,
                includes=includes,
                target_hint=extract_target_hint(output, directory),
            )
        )
    return commands


def parse_log(path: str) -> BuildSummary:
    summary = BuildSummary()
    current_context: Optional[str] = None
    current_tu: Optional[str] = None

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = strip_ansi(raw_line.rstrip("\n"))
            if not line.strip():
                continue

            if line.startswith("Change Dir: "):
                summary.build_dir = line.split(":", 1)[1].strip().strip("'")
                continue

            if line.startswith("Run Build Command(s): "):
                summary.build_command = line.split(":", 1)[1].strip()
                continue

            m = ENTER_DIR_RE.search(line)
            if m:
                summary.directories_entered.append(m.group(1))
                continue

            m = BUILDING_OBJ_RE.match(line)
            if m:
                summary.building_steps.append(m.group(1))
                continue

            m = LINKING_RE.match(line)
            if m:
                summary.linking_steps.append(m.group(1))
                continue

            m = BUILT_TARGET_RE.match(line)
            if m:
                summary.built_targets.append(m.group(1))
                continue

            cmd = parse_compile_command_from_line(line)
            if cmd:
                summary.compile_commands.append(cmd)
                current_tu = normalize_path(cmd.file)
                for wf in cmd.warning_flags:
                    summary.compiler_warning_flags_seen[wf] += 1
                if cmd.std:
                    summary.standards_seen[cmd.std] += 1
                if cmd.optimization:
                    summary.optimization_levels_seen[cmd.optimization] += 1
                continue

            m = IN_FUNCTION_RE.match(line)
            if m:
                current_context = compact_message(line.strip())
                continue

            m = INSTANTIATION_RE.match(line)
            if m:
                current_context = compact_message(line.strip())
                continue

            m = WARNING_RE.match(line)
            if m:
                d = Diagnostic(
                    kind="warning",
                    file=normalize_path(m.group("file")) or m.group("file"),
                    line=int(m.group("line")),
                    col=int(m.group("col")),
                    message=m.group("message").strip(),
                    flag=m.group("flag"),
                    context=current_context,
                    translation_unit=current_tu,
                )
                summary.diagnostics.append(d)
                summary.warnings += 1
                summary.warnings_by_file[d.file] += 1
                summary.warnings_by_flag[d.flag or "(unclassified)"] += 1
                continue

            m = ERROR_RE.match(line)
            if m:
                d = Diagnostic(
                    kind="error",
                    file=normalize_path(m.group("file")) or m.group("file"),
                    line=int(m.group("line")),
                    col=int(m.group("col")),
                    message=m.group("message").strip(),
                    flag=m.group("flag"),
                    context=current_context,
                    translation_unit=current_tu,
                )
                summary.diagnostics.append(d)
                summary.errors += 1
                summary.errors_by_file[d.file] += 1
                summary.raw_error_lines.append(line)
                continue

            m = NOTE_RE.match(line)
            if m:
                d = Diagnostic(
                    kind="note",
                    file=normalize_path(m.group("file")) or m.group("file"),
                    line=int(m.group("line")),
                    col=int(m.group("col")),
                    message=m.group("message").strip(),
                    context=current_context,
                    translation_unit=current_tu,
                )
                summary.diagnostics.append(d)
                summary.notes += 1
                continue

            if FATAL_RE.search(line):
                summary.fatal_lines.append(line)
            if MAKE_ERROR_RE.search(line) or FAILED_RE.search(line):
                summary.raw_error_lines.append(line)

    if summary.errors > 0 or summary.fatal_lines or summary.raw_error_lines:
        summary.success = False
    elif summary.built_targets:
        summary.success = True
    else:
        summary.success = None

    return summary


def merge_compile_database(summary: BuildSummary, compile_db: list[CompileCommand]) -> None:
    if not compile_db:
        return

    if not summary.compile_commands:
        summary.compile_commands = compile_db
    else:
        known = {(c.file, c.output, c.target_hint) for c in summary.compile_commands}
        for cmd in compile_db:
            key = (cmd.file, cmd.output, cmd.target_hint)
            if key not in known:
                summary.compile_commands.append(cmd)
                known.add(key)

    summary.compiler_warning_flags_seen.clear()
    summary.standards_seen.clear()
    summary.optimization_levels_seen.clear()

    for cmd in summary.compile_commands:
        for wf in cmd.warning_flags:
            summary.compiler_warning_flags_seen[wf] += 1
        if cmd.std:
            summary.standards_seen[cmd.std] += 1
        if cmd.optimization:
            summary.optimization_levels_seen[cmd.optimization] += 1

    tu_to_targets: dict[str, set[str]] = defaultdict(set)
    for cmd in summary.compile_commands:
        if cmd.file and cmd.target_hint:
            tu_to_targets[cmd.file].add(cmd.target_hint)

    for d in summary.diagnostics:
        if d.translation_unit and d.translation_unit in tu_to_targets:
            targets = sorted(tu_to_targets[d.translation_unit])
            if len(targets) == 1:
                d.target_hint = targets[0]
            elif targets:
                d.target_hint = ",".join(targets)


def unique_diagnostics(diags: list[Diagnostic], mode: str) -> list[Diagnostic]:
    seen = set()
    out = []
    for d in diags:
        key = d.site_key() if mode == "site" else d.broad_key()
        if key in seen:
            continue
        seen.add(key)
        out.append(d)
    return out


def top_counter(counter: Counter, n: int) -> list[tuple[str, int]]:
    return counter.most_common(n)


def unique_sources(summary: BuildSummary) -> int:
    return len({c.file for c in summary.compile_commands if c.file})


def format_status(summary: BuildSummary) -> str:
    if summary.success is True:
        return "SUCCESS"
    if summary.success is False:
        return "FAILED"
    return "UNKNOWN"


def is_header(path: str) -> bool:
    lower = path.lower()
    return lower.endswith((".h", ".hh", ".hpp", ".hxx"))


def diagnostic_repetition_map(diags: list[Diagnostic], mode: str) -> dict[tuple, int]:
    counter: Counter = Counter()
    for d in diags:
        key = d.site_key() if mode == "site" else d.broad_key()
        counter[key] += 1
    return dict(counter)


def summarize_warning_hotspots(summary: BuildSummary, dedupe_mode: str) -> dict:
    warnings = [d for d in summary.diagnostics if d.kind == "warning"]
    unique_warns = unique_diagnostics(warnings, dedupe_mode)
    raw_repeat_map = diagnostic_repetition_map(warnings, dedupe_mode)

    unique_by_file = Counter()
    raw_by_file = Counter()
    unique_by_flag = Counter()
    raw_by_flag = Counter()
    repeated_sites = []

    for d in warnings:
        raw_by_file[d.file] += 1
        raw_by_flag[d.flag or "(unclassified)"] += 1

    for d in unique_warns:
        unique_by_file[d.file] += 1
        unique_by_flag[d.flag or "(unclassified)"] += 1
        key = d.site_key() if dedupe_mode == "site" else d.broad_key()
        repeats = raw_repeat_map.get(key, 1)
        if repeats > 1:
            repeated_sites.append((d, repeats))

    repeated_sites.sort(key=lambda item: item[1], reverse=True)

    header_raw = sum(1 for d in warnings if is_header(d.file))
    header_unique = sum(1 for d in unique_warns if is_header(d.file))
    source_raw = len(warnings) - header_raw
    source_unique = len(unique_warns) - header_unique

    return {
        "raw_count": len(warnings),
        "unique_count": len(unique_warns),
        "raw_by_file": raw_by_file,
        "unique_by_file": unique_by_file,
        "raw_by_flag": raw_by_flag,
        "unique_by_flag": unique_by_flag,
        "header_raw": header_raw,
        "header_unique": header_unique,
        "source_raw": source_raw,
        "source_unique": source_unique,
        "repeated_sites": repeated_sites,
        "unique_warnings": unique_warns,
    }


def representative_diagnostics(diags: list[Diagnostic], max_examples: int, dedupe_mode: str) -> list[Diagnostic]:
    return unique_diagnostics(diags, dedupe_mode)[:max_examples]


def format_diag_text_lines(d: Diagnostic, build_dir: Optional[str]) -> list[str]:
    location = f"{shorten_path(d.file, build_dir)}:{d.line}:{d.col}"
    flag = f" [{d.flag}]" if d.flag else ""
    lines = [f"- {location}{flag}", f"  {compact_message(d.message)}"]
    if d.context:
        lines.append(f"  Context: {compact_message(d.context)}")
    if d.translation_unit:
        lines.append(f"  TU: {shorten_path(d.translation_unit, build_dir)}")
    if d.target_hint:
        lines.append(f"  Target: {d.target_hint}")
    return lines


def format_diag_markdown_lines(d: Diagnostic, build_dir: Optional[str]) -> list[str]:
    location = f"{shorten_path(d.file, build_dir)}:{d.line}:{d.col}"
    flag = f" `{d.flag}`" if d.flag else ""
    lines = [f"- **{location}**{flag}: {compact_message(d.message)}"]
    if d.context:
        lines.append(f"  - Context: `{compact_message(d.context)}`")
    if d.translation_unit:
        lines.append(f"  - TU: `{shorten_path(d.translation_unit, build_dir)}`")
    if d.target_hint:
        lines.append(f"  - Target: `{d.target_hint}`")
    return lines


def generate_human_summary(
    summary: BuildSummary,
    dedupe_mode: str = "site",
    max_examples: int = 5,
    show_files: int = 10,
    show_flags: int = 10,
    show_repeated: int = 5,
) -> str:
    lines: list[str] = []
    warn_stats = summarize_warning_hotspots(summary, dedupe_mode)
    warnings = [d for d in summary.diagnostics if d.kind == "warning"]
    errors = [d for d in summary.diagnostics if d.kind == "error"]

    lines.append("BUILD SUMMARY")
    lines.append("=" * 72)
    lines.append(f"Status           : {format_status(summary)}")
    if summary.build_dir:
        lines.append(f"Build directory  : {shorten_path(summary.build_dir, summary.build_dir)}")
    if summary.build_command:
        lines.append(f"Build command    : {summary.build_command}")
    lines.append(f"Dedupe mode      : {dedupe_mode}")
    lines.append("")

    lines.append("High-level results")
    lines.append("-" * 72)
    lines.append(f"Raw warnings     : {warn_stats['raw_count']}")
    lines.append(f"Unique warnings  : {warn_stats['unique_count']}")
    lines.append(f"Errors           : {summary.errors}")
    lines.append(f"Notes            : {summary.notes}")
    lines.append(f"Compile commands : {len(summary.compile_commands)}")
    lines.append(f"Unique sources   : {unique_sources(summary)}")
    lines.append(f"Built targets    : {len(summary.built_targets)}")
    lines.append("")

    lines.append("Warning origin breakdown")
    lines.append("-" * 72)
    lines.append(f"Header warnings  : raw={warn_stats['header_raw']}, unique={warn_stats['header_unique']}")
    lines.append(f"Source warnings  : raw={warn_stats['source_raw']}, unique={warn_stats['source_unique']}")
    if warn_stats["raw_count"] > 0:
        lines.append(f"Repeat factor    : {warn_stats['raw_count'] / max(warn_stats['unique_count'], 1):.2f}x raw-to-unique")
    lines.append("")

    if summary.built_targets:
        lines.append("Targets built")
        lines.append("-" * 72)
        for t in summary.built_targets[:20]:
            lines.append(f"- {t}")
        if len(summary.built_targets) > 20:
            lines.append(f"... and {len(summary.built_targets) - 20} more")
        lines.append("")

    if summary.standards_seen or summary.optimization_levels_seen or summary.compiler_warning_flags_seen:
        lines.append("Compiler configuration detected")
        lines.append("-" * 72)
        if summary.standards_seen:
            lines.append("C++ standards    : " + ", ".join(f"{k} ({v})" for k, v in summary.standards_seen.most_common()))
        if summary.optimization_levels_seen:
            lines.append("Optimization     : " + ", ".join(f"{k} ({v})" for k, v in summary.optimization_levels_seen.most_common()))
        if summary.compiler_warning_flags_seen:
            top_flags_text = ", ".join(
                f"{flag} ({count})"
                for flag, count in summary.compiler_warning_flags_seen.most_common(min(12, len(summary.compiler_warning_flags_seen)))
            )
            lines.append(f"Warning flags    : {top_flags_text}")
        lines.append("")

    lines.append("Most common warning categories")
    lines.append("-" * 72)
    for flag, count in top_counter(warn_stats["raw_by_flag"], show_flags):
        uniq = warn_stats["unique_by_flag"].get(flag, 0)
        lines.append(f"- {flag}: raw={count}, unique={uniq}")
    lines.append("")

    lines.append("Files with the most warnings")
    lines.append("-" * 72)
    for file, count in top_counter(warn_stats["raw_by_file"], show_files):
        uniq = warn_stats["unique_by_file"].get(file, 0)
        lines.append(f"- {shorten_path(file, summary.build_dir)}: raw={count}, unique={uniq}")
    lines.append("")

    if warn_stats["repeated_sites"]:
        lines.append("Most repeated warning sites")
        lines.append("-" * 72)
        for d, repeats in warn_stats["repeated_sites"][:show_repeated]:
            location = f"{shorten_path(d.file, summary.build_dir)}:{d.line}:{d.col}"
            flag = f" [{d.flag}]" if d.flag else ""
            lines.append(f"- {location}{flag}: repeated {repeats} times")
            lines.append(f"  {compact_message(d.message)}")
        lines.append("")

    warn_examples = representative_diagnostics(warnings, max_examples, dedupe_mode)
    if warn_examples:
        lines.append("Representative warnings")
        lines.append("-" * 72)
        for d in warn_examples:
            lines.extend(format_diag_text_lines(d, summary.build_dir))
        lines.append("")

    err_examples = representative_diagnostics(errors, max_examples, dedupe_mode)
    if err_examples:
        lines.append("Representative errors")
        lines.append("-" * 72)
        for d in err_examples:
            lines.extend(format_diag_text_lines(d, summary.build_dir))
        lines.append("")

    lines.append("Interpretation")
    lines.append("-" * 72)
    if summary.success is True and summary.errors == 0:
        if warn_stats["raw_count"] == 0:
            lines.append("The build completed cleanly with no warnings or errors.")
        else:
            top_raw = warn_stats["raw_by_flag"].most_common(3)
            top_text = ", ".join(f"{flag} ({count})" for flag, count in top_raw)
            lines.append(f"The build succeeded, but warnings remain. Dominant raw warning classes: {top_text}.")
            if warn_stats["raw_count"] != warn_stats["unique_count"]:
                lines.append(
                    f"Many warnings are repeated emissions: {warn_stats['raw_count']} raw warnings collapse to "
                    f"{warn_stats['unique_count']} unique warning sites."
                )
            if warn_stats["header_raw"] > warn_stats["source_raw"]:
                lines.append(
                    "Most warning volume originates in headers, which suggests a small number of warning sites "
                    "are being re-emitted across many translation units."
                )
    elif summary.success is False:
        if summary.errors > 0:
            lines.append("The build failed because compiler errors were detected.")
        else:
            lines.append("The build appears to have failed, but no structured compiler error lines were parsed.")
    else:
        lines.append("The parser could not conclusively determine build status.")

    top_flag = warn_stats["raw_by_flag"].most_common(1)
    if top_flag:
        flag_name, _ = top_flag[0]
        if flag_name == "-Wshadow":
            lines.append("The biggest cleanup opportunity is shadowed identifiers.")
        elif flag_name in {"-Wconversion", "-Wsign-conversion"}:
            lines.append("The biggest cleanup opportunity is numeric conversion hygiene.")

    return "\n".join(lines)


def generate_markdown_summary(
    summary: BuildSummary,
    dedupe_mode: str = "site",
    max_examples: int = 5,
    show_files: int = 10,
    show_flags: int = 10,
    show_repeated: int = 5,
) -> str:
    warn_stats = summarize_warning_hotspots(summary, dedupe_mode)
    warnings = [d for d in summary.diagnostics if d.kind == "warning"]
    errors = [d for d in summary.diagnostics if d.kind == "error"]

    lines: list[str] = []
    lines.append("# Build Summary")
    lines.append("")
    lines.append("| Field | Value |")
    lines.append("|---|---|")
    lines.append(f"| Status | {format_status(summary)} |")
    if summary.build_dir:
        lines.append(f"| Build directory | `{shorten_path(summary.build_dir, summary.build_dir)}` |")
    if summary.build_command:
        lines.append(f"| Build command | `{summary.build_command}` |")
    lines.append(f"| Dedupe mode | `{dedupe_mode}` |")
    lines.append("")

    lines.append("## High-level Results")
    lines.append("")
    lines.append("| Metric | Value |")
    lines.append("|---|---:|")
    lines.append(f"| Raw warnings | {warn_stats['raw_count']} |")
    lines.append(f"| Unique warnings | {warn_stats['unique_count']} |")
    lines.append(f"| Errors | {summary.errors} |")
    lines.append(f"| Notes | {summary.notes} |")
    lines.append(f"| Compile commands | {len(summary.compile_commands)} |")
    lines.append(f"| Unique sources | {unique_sources(summary)} |")
    lines.append(f"| Built targets | {len(summary.built_targets)} |")
    lines.append("")

    lines.append("## Warning Origin Breakdown")
    lines.append("")
    lines.append("| Category | Raw | Unique |")
    lines.append("|---|---:|---:|")
    lines.append(f"| Headers | {warn_stats['header_raw']} | {warn_stats['header_unique']} |")
    lines.append(f"| Sources | {warn_stats['source_raw']} | {warn_stats['source_unique']} |")
    if warn_stats["raw_count"] > 0:
        lines.append("")
        lines.append(f"- Repeat factor: **{warn_stats['raw_count'] / max(warn_stats['unique_count'], 1):.2f}x** raw-to-unique")
    lines.append("")

    if summary.built_targets:
        lines.append("## Targets Built")
        lines.append("")
        for t in summary.built_targets[:20]:
            lines.append(f"- `{t}`")
        if len(summary.built_targets) > 20:
            lines.append(f"- ... and {len(summary.built_targets) - 20} more")
        lines.append("")

    if summary.standards_seen or summary.optimization_levels_seen or summary.compiler_warning_flags_seen:
        lines.append("## Compiler Configuration")
        lines.append("")
        if summary.standards_seen:
            lines.append("- C++ standards: " + ", ".join(f"`{k}` ({v})" for k, v in summary.standards_seen.most_common()))
        if summary.optimization_levels_seen:
            lines.append("- Optimization: " + ", ".join(f"`{k}` ({v})" for k, v in summary.optimization_levels_seen.most_common()))
        if summary.compiler_warning_flags_seen:
            lines.append(
                "- Warning flags: "
                + ", ".join(
                    f"`{flag}` ({count})"
                    for flag, count in summary.compiler_warning_flags_seen.most_common(min(12, len(summary.compiler_warning_flags_seen)))
                )
            )
        lines.append("")

    lines.append("## Most Common Warning Categories")
    lines.append("")
    lines.append("| Warning flag | Raw | Unique |")
    lines.append("|---|---:|---:|")
    for flag, count in top_counter(warn_stats["raw_by_flag"], show_flags):
        uniq = warn_stats["unique_by_flag"].get(flag, 0)
        lines.append(f"| `{flag}` | {count} | {uniq} |")
    lines.append("")

    lines.append("## Files with the Most Warnings")
    lines.append("")
    lines.append("| File | Raw | Unique |")
    lines.append("|---|---:|---:|")
    for file, count in top_counter(warn_stats["raw_by_file"], show_files):
        uniq = warn_stats["unique_by_file"].get(file, 0)
        lines.append(f"| `{shorten_path(file, summary.build_dir)}` | {count} | {uniq} |")
    lines.append("")

    if warn_stats["repeated_sites"]:
        lines.append("## Most Repeated Warning Sites")
        lines.append("")
        for d, repeats in warn_stats["repeated_sites"][:show_repeated]:
            location = f"{shorten_path(d.file, summary.build_dir)}:{d.line}:{d.col}"
            flag = f" `{d.flag}`" if d.flag else ""
            lines.append(f"- **{location}**{flag}: repeated **{repeats}** times")
            lines.append(f"  - {compact_message(d.message)}")
        lines.append("")

    warn_examples = representative_diagnostics(warnings, max_examples, dedupe_mode)
    if warn_examples:
        lines.append("## Representative Warnings")
        lines.append("")
        for d in warn_examples:
            lines.extend(format_diag_markdown_lines(d, summary.build_dir))
        lines.append("")

    err_examples = representative_diagnostics(errors, max_examples, dedupe_mode)
    if err_examples:
        lines.append("## Representative Errors")
        lines.append("")
        for d in err_examples:
            lines.extend(format_diag_markdown_lines(d, summary.build_dir))
        lines.append("")

    lines.append("## Interpretation")
    lines.append("")
    if summary.success is True and summary.errors == 0:
        if warn_stats["raw_count"] == 0:
            lines.append("- The build completed cleanly with no warnings or errors.")
        else:
            top_raw = warn_stats["raw_by_flag"].most_common(3)
            top_text = ", ".join(f"`{flag}` ({count})" for flag, count in top_raw)
            lines.append(f"- The build succeeded, but warnings remain. Dominant raw warning classes: {top_text}.")
            if warn_stats["raw_count"] != warn_stats["unique_count"]:
                lines.append(
                    f"- Many warnings are repeated emissions: **{warn_stats['raw_count']}** raw warnings collapse to "
                    f"**{warn_stats['unique_count']}** unique warning sites."
                )
            if warn_stats["header_raw"] > warn_stats["source_raw"]:
                lines.append(
                    "- Most warning volume originates in headers, which suggests a small number of warning sites are "
                    "being re-emitted across many translation units."
                )
    elif summary.success is False:
        if summary.errors > 0:
            lines.append("- The build failed because compiler errors were detected.")
        else:
            lines.append("- The build appears to have failed, but no structured compiler error lines were parsed.")
    else:
        lines.append("- The parser could not conclusively determine build status.")

    top_flag = warn_stats["raw_by_flag"].most_common(1)
    if top_flag:
        flag_name, _ = top_flag[0]
        if flag_name == "-Wshadow":
            lines.append("- The biggest cleanup opportunity is shadowed identifiers.")
        elif flag_name in {"-Wconversion", "-Wsign-conversion"}:
            lines.append("- The biggest cleanup opportunity is numeric conversion hygiene.")

    return "\n".join(lines)


def summary_to_json(summary: BuildSummary, dedupe_mode: str, max_examples: int) -> dict:
    warn_stats = summarize_warning_hotspots(summary, dedupe_mode)
    warnings = [d for d in summary.diagnostics if d.kind == "warning"]
    errors = [d for d in summary.diagnostics if d.kind == "error"]

    return {
        "status": format_status(summary),
        "build_dir": summary.build_dir,
        "build_command": summary.build_command,
        "dedupe_mode": dedupe_mode,
        "counts": {
            "warnings_raw": warn_stats["raw_count"],
            "warnings_unique": warn_stats["unique_count"],
            "errors": summary.errors,
            "notes": summary.notes,
            "compile_commands": len(summary.compile_commands),
            "unique_sources": unique_sources(summary),
            "built_targets": len(summary.built_targets),
            "header_warnings_raw": warn_stats["header_raw"],
            "header_warnings_unique": warn_stats["header_unique"],
            "source_warnings_raw": warn_stats["source_raw"],
            "source_warnings_unique": warn_stats["source_unique"],
        },
        "built_targets": summary.built_targets,
        "warning_categories_raw": dict(warn_stats["raw_by_flag"].most_common()),
        "warning_categories_unique": dict(warn_stats["unique_by_flag"].most_common()),
        "warnings_by_file_raw": dict(warn_stats["raw_by_file"].most_common()),
        "warnings_by_file_unique": dict(warn_stats["unique_by_file"].most_common()),
        "errors_by_file": dict(summary.errors_by_file.most_common()),
        "compiler_config": {
            "standards_seen": dict(summary.standards_seen.most_common()),
            "optimization_levels_seen": dict(summary.optimization_levels_seen.most_common()),
            "warning_flags_seen_in_commands": dict(summary.compiler_warning_flags_seen.most_common()),
        },
        "most_repeated_warning_sites": [
            {"repeats": repeats, "diagnostic": asdict(d)}
            for d, repeats in warn_stats["repeated_sites"][:max_examples]
        ],
        "examples": {
            "warnings": [asdict(d) for d in representative_diagnostics(warnings, max_examples, dedupe_mode)],
            "errors": [asdict(d) for d in representative_diagnostics(errors, max_examples, dedupe_mode)],
        },
    }


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Summarize a C/C++ build log.")
    p.add_argument("logfile", help="Path to build log")
    p.add_argument("--compile-db", help="Optional path to compile_commands.json")
    p.add_argument("--json", action="store_true", help="Emit JSON output")
    p.add_argument("--markdown", action="store_true", help="Emit Markdown output")
    p.add_argument("--max-examples", type=int, default=5)
    p.add_argument("--show-files", type=int, default=10)
    p.add_argument("--show-flags", type=int, default=10)
    p.add_argument("--show-repeated", type=int, default=5)
    p.add_argument(
        "--dedupe-mode",
        choices=["site", "broad"],
        default="site",
        help="site=(file,line,col,message,flag), broad=(file,line,message,flag)",
    )
    return p


def main() -> int:
    args = build_arg_parser().parse_args()

    if args.json and args.markdown:
        print("error: choose only one of --json or --markdown", file=sys.stderr)
        return 2

    if not os.path.isfile(args.logfile):
        print(f"error: file not found: {args.logfile}", file=sys.stderr)
        return 2

    summary = parse_log(args.logfile)

    if args.compile_db:
        if not os.path.isfile(args.compile_db):
            print(f"error: compile database not found: {args.compile_db}", file=sys.stderr)
            return 2
        compile_db = load_compile_database(args.compile_db)
        merge_compile_database(summary, compile_db)

    if args.json:
        print(json.dumps(summary_to_json(summary, args.dedupe_mode, args.max_examples), indent=2))
    elif args.markdown:
        print(
            generate_markdown_summary(
                summary,
                dedupe_mode=args.dedupe_mode,
                max_examples=args.max_examples,
                show_files=args.show_files,
                show_flags=args.show_flags,
                show_repeated=args.show_repeated,
            )
        )
    else:
        print(
            generate_human_summary(
                summary,
                dedupe_mode=args.dedupe_mode,
                max_examples=args.max_examples,
                show_files=args.show_files,
                show_flags=args.show_flags,
                show_repeated=args.show_repeated,
            )
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
