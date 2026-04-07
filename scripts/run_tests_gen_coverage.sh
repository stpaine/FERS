#!/bin/bash

set -Eeuo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly SCRIPT_DIR
readonly REPO_ROOT="${SCRIPT_DIR}/.."
readonly COVERAGE_DIR="${REPO_ROOT}/build/coverage"

die() {
	printf 'Error: %s\n' "$*" >&2
	exit 1
}

on_err() {
	local exit_code=$?
	printf 'Error: command failed (exit %d) at %s:%s: %s\n' \
		"$exit_code" "${BASH_SOURCE[0]}" "${BASH_LINENO[0]}" "${BASH_COMMAND}" >&2
	exit "$exit_code"
}

require_cmd() {
	command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

trap on_err ERR

for cmd in cmake ctest lcov genhtml; do
	require_cmd "$cmd"
done

cd "$REPO_ROOT"

cmake --preset=coverage
cmake --build --preset=coverage --parallel
ctest --preset=coverage --parallel

cd "$COVERAGE_DIR"

lcov --capture --directory . --output-file coverage.info --ignore-errors source,gcov,negative
lcov --remove coverage.info '/usr/*' '*/tests/*' '*/vcpkg_installed/*' --output-file coverage_filtered.info
genhtml coverage_filtered.info --output-directory coverage_report --ignore-errors source
