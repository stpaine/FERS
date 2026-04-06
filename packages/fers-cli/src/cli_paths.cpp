// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

#include "cli_paths.h"

namespace core
{
	std::filesystem::path resolveOutputDir(const std::string& script_file,
										   const std::optional<std::string>& output_dir) noexcept
	{
		if (output_dir)
		{
			return std::filesystem::path(*output_dir);
		}

		std::filesystem::path default_out_dir = std::filesystem::path(script_file).parent_path();
		if (default_out_dir.empty())
		{
			default_out_dir = ".";
		}

		return default_out_dir;
	}

	std::filesystem::path resolveKmlOutputPath(const std::string& script_file,
											   const std::filesystem::path& final_output_dir,
											   const std::optional<std::string>& kml_file) noexcept
	{
		if (kml_file && !kml_file->empty())
		{
			const std::filesystem::path provided_kml_path(*kml_file);
			if (provided_kml_path.has_parent_path() || provided_kml_path.is_absolute())
			{
				return provided_kml_path;
			}

			return final_output_dir / provided_kml_path;
		}

		std::filesystem::path kml_output_path = final_output_dir / std::filesystem::path(script_file).filename();
		kml_output_path.replace_extension(".kml");
		return kml_output_path;
	}
}
