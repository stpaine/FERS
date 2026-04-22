#!/usr/bin/env python3
"""Verify that repo-managed version fields stay in sync."""

from __future__ import annotations

import json
import re
import sys
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_toml(path: Path) -> dict:
    return tomllib.loads(path.read_text(encoding="utf-8"))


def read_bun_workspace_version(path: Path, workspace_path: str) -> str | None:
    lock_text = path.read_text(encoding="utf-8")
    workspace_match = re.search(
        rf'"{re.escape(workspace_path)}":\s*\{{(?P<body>.*?)^\s*\}}',
        lock_text,
        re.MULTILINE | re.DOTALL,
    )
    if workspace_match is None:
        return None

    version_match = re.search(r'^\s*"version":\s*"([^"]+)"', workspace_match.group("body"), re.MULTILINE)
    if version_match is None:
        return None

    return version_match.group(1)


def main() -> int:
    expected = (ROOT / "version.txt").read_text(encoding="utf-8").strip()
    errors: list[str] = []

    root_package = read_json(ROOT / "package.json")
    if root_package.get("version") != expected:
        errors.append("package.json version does not match version.txt")

    vcpkg = read_json(ROOT / "vcpkg.json")
    if vcpkg.get("version-semver") != expected:
        errors.append("vcpkg.json version-semver does not match version.txt")

    ui_package = read_json(ROOT / "packages/fers-ui/package.json")
    if ui_package.get("version") != expected:
        errors.append("packages/fers-ui/package.json version does not match version.txt")

    bun_workspace_version = read_bun_workspace_version(ROOT / "bun.lock", "packages/fers-ui")
    if bun_workspace_version != expected:
        errors.append("bun.lock fers-ui workspace version does not match version.txt")

    tauri_cargo = read_toml(ROOT / "packages/fers-ui/src-tauri/Cargo.toml")
    if tauri_cargo.get("package", {}).get("version") != expected:
        errors.append("packages/fers-ui/src-tauri/Cargo.toml package.version does not match version.txt")

    tauri_config = read_json(ROOT / "packages/fers-ui/src-tauri/tauri.conf.json")
    if tauri_config.get("version") != "../package.json":
        errors.append("packages/fers-ui/src-tauri/tauri.conf.json version must point to ../package.json")

    root_cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    if "version.txt" not in root_cmake or "FERS_PROJECT_VERSION" not in root_cmake:
        errors.append("CMakeLists.txt is not wired to read the shared version from version.txt")

    if errors:
        print("Version verification failed:", file=sys.stderr)
        for error in errors:
            print(f" - {error}", file=sys.stderr)
        return 1

    print(f"All managed versions match {expected}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
