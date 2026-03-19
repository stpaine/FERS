// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file hdf5_handler.h
 * @brief Header file for HDF5 data export and import functions.
 */

#pragma once

#include <string>
#include <vector>

#include "core/config.h"

namespace HighFive
{
	class File;
}

namespace serial
{
	/**
	 * @brief Adds a chunk of data to an HDF5 file.
	 *
	 * @param file The HDF5 file where the chunk is written.
	 * @param data A vector of complex data to be written.
	 * @param time The time attribute associated with the chunk.
	 * @param fullscale The fullscale attribute for the chunk.
	 * @param count The sequential count number for chunk naming.
	 * @throws std::runtime_error If there is an error writing data or setting attributes.
	 */
	void addChunkToFile(HighFive::File& file, const std::vector<ComplexType>& data, RealType time, RealType fullscale,
						unsigned count);

	/**
	 * @brief Reads pulse data from an HDF5 file.
	 *
	 * @param name The name of the HDF5 file.
	 * @param data A reference to a vector where the complex data will be stored.
	 * @throws std::runtime_error If the file does not exist or the datasets "I" and "Q" have mismatched sizes.
	 */
	void readPulseData(const std::string& name, std::vector<ComplexType>& data);

	/**
	 * @brief Reads a 2D pattern dataset from an HDF5 file.
	 *
	 * @param name The name of the HDF5 file.
	 * @param datasetName The name of the dataset to be read.
	 * @return A 2D vector containing the pattern data.
	 * @throws std::runtime_error If there is an error handling the file or if the dataset dimensions are invalid.
	 */
	std::vector<std::vector<RealType>> readPattern(const std::string& name, const std::string& datasetName);
}
