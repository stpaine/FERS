// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file waveform_factory.cpp
 * @brief Implementation for loading waveform data into RadarSignal objects.
 */

#include "waveform_factory.h"

#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "core/config.h"
#include "core/parameters.h"
#include "core/sim_id.h"
#include "hdf5_handler.h"
#include "signal/radar_signal.h"

using fers_signal::RadarSignal;
using fers_signal::Signal;

namespace
{
	/**
	 * @brief Loads a radar waveform from an HDF5 file and returns a RadarSignal object.
	 *
	 * @param name The name of the radar signal.
	 * @param filepath The path to the HDF5 file containing the waveform data.
	 * @param power The power of the radar signal in the waveform.
	 * @param carrierFreq The carrier frequency of the radar signal.
	 * @return A unique pointer to a RadarSignal object loaded with the waveform data.
	 * @throws std::runtime_error If the file cannot be opened or the file format is unrecognized.
	 */
	std::unique_ptr<RadarSignal> loadWaveformFromHdf5File(const std::string& name,
														  const std::filesystem::path& filepath, const RealType power,
														  const RealType carrierFreq, const SimId id)
	{
		std::vector<ComplexType> data;
		serial::readPulseData(filepath.string(), data);

		auto signal = std::make_unique<Signal>();
		signal->load(data, data.size(), params::rate());
		return std::make_unique<RadarSignal>(
			name, power, carrierFreq, static_cast<RealType>(data.size()) / params::rate(), std::move(signal), id);
	}

	/**
	 * @brief Loads a radar waveform from a CSV file and returns a RadarSignal object.
	 *
	 * @param name The name of the radar signal.
	 * @param filepath The path to the CSV file containing the waveform data.
	 * @param power The power of the radar signal in the waveform.
	 * @param carrierFreq The carrier frequency of the radar signal.
	 * @return A unique pointer to a RadarSignal object loaded with the waveform data.
	 * @throws std::runtime_error If the file cannot be opened or the file format is unrecognized.
	 */
	std::unique_ptr<RadarSignal> loadWaveformFromCsvFile(const std::string& name, const std::filesystem::path& filepath,
														 const RealType power, const RealType carrierFreq,
														 const SimId id)
	{
		std::ifstream ifile(filepath);
		if (!ifile)
		{
			LOG(logging::Level::FATAL, "Could not open file '{}' to read waveform", filepath.string());
			throw std::runtime_error("Could not open file '" + filepath.string() + "' to read waveform");
		}

		RealType rlength, rate;
		ifile >> rlength >> rate;

		const auto length = static_cast<std::size_t>(rlength);
		std::vector<ComplexType> data(length);

		// Read the file data
		for (std::size_t done = 0; done < length && ifile >> data[done]; ++done)
		{
		}

		if (ifile.fail() || data.size() != length)
		{
			LOG(logging::Level::FATAL, "Could not read full waveform from file '{}'", filepath.string());
			throw std::runtime_error("Could not read full waveform from file '" + filepath.string() + "'");
		}

		auto signal = std::make_unique<Signal>();
		signal->load(data, length, rate);
		return std::make_unique<RadarSignal>(name, power, carrierFreq, rlength / rate, std::move(signal), id);
	}

	/**
	 * @brief Checks if a filename has a specific extension.
	 *
	 * @param filename The filename to check.
	 * @param ext The extension to check for.
	 * @return True if the filename has the specified extension, false otherwise.
	 */
	constexpr bool hasExtension(const std::string_view filename, const std::string_view ext) noexcept
	{
		return filename.ends_with(ext);
	}
}

namespace serial
{
	std::unique_ptr<RadarSignal> loadWaveformFromFile(const std::string& name, const std::string& filename,
													  const RealType power, const RealType carrierFreq, const SimId id)
	{
		const std::filesystem::path filepath = filename;
		const auto extension = filepath.extension().string();

		if (hasExtension(extension, ".csv"))
		{
			auto wave = loadWaveformFromCsvFile(name, filepath, power, carrierFreq, id);
			wave->setFilename(filename);
			return wave;
		}
		if (hasExtension(extension, ".h5"))
		{
			auto wave = loadWaveformFromHdf5File(name, filepath, power, carrierFreq, id);
			wave->setFilename(filename);
			return wave;
		}

		LOG(logging::Level::FATAL, "Unrecognized file extension '{}' for file: '{}'", extension, filename);
		throw std::runtime_error("Unrecognized file extension '" + extension + "' for file: " + filename);
	}
}
