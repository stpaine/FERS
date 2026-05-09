// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

/**
 * @file main.cpp
 * @brief Entry point for the FERS command-line interface (CLI).
 *
 * This executable acts as a wrapper around the libfers core library. It parses
 * command-line arguments, uses the libfers C-API to load and run a simulation,
 * and reports progress to the console.
 */

#include "cli_runner.h"

int main(const int argc, char* argv[]) { return core::runCli(argc, argv); }
