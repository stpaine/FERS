#!/bin/bash

set -Eeuo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly SCRIPT_DIR
readonly REPO_ROOT="${SCRIPT_DIR}/.."
readonly COVERAGE_DIR="${REPO_ROOT}/build/coverage"
# lcov suppresses a warning category only when the category appears twice.
readonly LCOV_IGNORE_ERRORS="source,gcov,negative,inconsistent,inconsistent"
readonly GENHTML_IGNORE_ERRORS="source,inconsistent,inconsistent"

TESTS_ONLY=false
SHOW_PASSING=false

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

print_failed_ctest_output() {
	local log_file="$1"

	awk '
		function is_result_line(line) {
			return line ~ /^[[:space:]]*[0-9]+\/[0-9]+ Test[[:space:]]+#[0-9]+:/
		}

		function is_boundary(line) {
			return line ~ /^Test project / ||
				line ~ /^[[:space:]]*Start[[:space:]][0-9]+:/ ||
				is_result_line(line) ||
				line ~ /^[[:space:]]*[0-9]+% tests passed/ ||
				line ~ /^Total Test time/ ||
				line ~ /^The following tests FAILED:/ ||
				line ~ /^Errors while running CTest/
		}

		/\*\*\*(Failed|Timeout|Not Run|Exception|SEGFAULT|Subprocess aborted)/ {
			if (printed) {
				print ""
			}
			print
			capturing_failure = 1
			printed = 1
			next
		}

		/^The following tests FAILED:/ {
			if (printed) {
				print ""
			}
			print
			capturing_failure = 0
			capturing_failed_list = 1
			printed = 1
			next
		}

		capturing_failed_list {
			if ($0 ~ /^Errors while running CTest/) {
				next
			}
			print
			next
		}

		capturing_failure {
			if (is_boundary($0)) {
				capturing_failure = 0
				next
			}
			print
		}

		END {
			exit printed ? 0 : 1
		}
	' "$log_file" >&2
}

run_ctest_failed_only() {
	local log_file
	local exit_code

	log_file="$(mktemp "${TMPDIR:-/tmp}/fers-ctest.XXXXXX")" || die "failed to create CTest log file"

	set +e
	ctest --preset=coverage --parallel --output-on-failure >"$log_file" 2>&1
	exit_code=$?
	set -e

	if [ "$exit_code" -ne 0 ]; then
		print_failed_ctest_output "$log_file" || {
			printf 'CTest failed before reporting failed tests:\n' >&2
			sed -n '1,200p' "$log_file" >&2
		}
	fi

	rm -f "$log_file"
	return "$exit_code"
}

run_ctest_show_passing() {
	local exit_code

	set +e
	ctest --preset=coverage --parallel --output-on-failure
	exit_code=$?
	set -e

	return "$exit_code"
}

run_ctest() {
	if [ "$SHOW_PASSING" = true ]; then
		run_ctest_show_passing
	else
		run_ctest_failed_only
	fi
}

# --- argument parsing ---
for arg in "$@"; do
	case "$arg" in
		--tests-only)
			TESTS_ONLY=true
			;;
		--show-passing)
			SHOW_PASSING=true
			;;
		*)
			die "unknown argument: $arg"
			;;
	esac
done

trap on_err ERR

# Only require lcov/genhtml if we actually use them
for cmd in cmake ctest; do
	require_cmd "$cmd"
done

if [ "$SHOW_PASSING" = false ]; then
	for cmd in awk mktemp sed; do
		require_cmd "$cmd"
	done
fi

if [ "$TESTS_ONLY" = false ]; then
	for cmd in lcov genhtml; do
		require_cmd "$cmd"
	done
fi

cd "$REPO_ROOT"

cmake --preset=coverage
cmake --build --preset=coverage --parallel
run_ctest || exit "$?"

# --- early exit if tests-only ---
if [ "$TESTS_ONLY" = true ]; then
	exit 0
fi

cd "$COVERAGE_DIR"

lcov --capture --directory . --output-file coverage.info \
	--rc geninfo_unexecuted_blocks=1 \
	--ignore-errors "$LCOV_IGNORE_ERRORS"
lcov --remove coverage.info '/usr/*' '*/tests/*' '*/vcpkg_installed/*' --output-file coverage_filtered.info \
	--ignore-errors "$LCOV_IGNORE_ERRORS"
genhtml coverage_filtered.info --output-directory coverage_report \
	--ignore-errors "$GENHTML_IGNORE_ERRORS"
