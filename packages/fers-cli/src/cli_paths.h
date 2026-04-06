// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

#pragma once

#include <filesystem>
#include <optional>
#include <string>

namespace core
{
	std::filesystem::path resolveOutputDir(const std::string& script_file,
										   const std::optional<std::string>& output_dir) noexcept;

	std::filesystem::path resolveKmlOutputPath(const std::string& script_file,
											   const std::filesystem::path& final_output_dir,
											   const std::optional<std::string>& kml_file) noexcept;
}
