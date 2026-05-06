// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

#pragma once

#include <filesystem>
#include <optional>
#include <string>

namespace core
{
	/// Resolves the final simulation output directory from CLI arguments.
	std::filesystem::path resolveOutputDir(const std::string& script_file,
										   const std::optional<std::string>& output_dir) noexcept;

	/// Resolves the KML output file path from CLI arguments and output directory.
	std::filesystem::path resolveKmlOutputPath(const std::string& script_file,
											   const std::filesystem::path& final_output_dir,
											   const std::optional<std::string>& kml_file) noexcept;
}
