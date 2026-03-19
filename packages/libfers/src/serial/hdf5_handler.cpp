// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file hdf5_handler.cpp
 * @brief Source file for HDF5 data export and import functions.
 */

#include "hdf5_handler.h"

#include <algorithm>
#include <complex>
#include <filesystem>
#include <format>
#include <highfive/highfive.hpp>
#include <stdexcept>

#include "core/logging.h"
#include "core/parameters.h"

using logging::Level;

namespace serial
{
	void readPulseData(const std::string& name, std::vector<ComplexType>& data)
	{
		if (!std::filesystem::exists(name))
		{
			LOG(Level::FATAL, "File '{}' not found", name);
			throw std::runtime_error("File " + name + " not found.");
		}

		LOG(Level::TRACE, "Opening file '{}'", name);
		const HighFive::File file(name, HighFive::File::ReadOnly);

		// Helper lambda to open group and read dataset
		auto read_dataset = [&file](const std::string& groupName, std::vector<double>& buffer) -> size_t
		{
			const auto group = file.getGroup("/" + groupName);

			const auto dataset = group.getDataSet("value");

			const auto dimensions = dataset.getSpace().getDimensions();
			const auto size = dimensions[0];

			buffer.resize(size);
			dataset.read(buffer);

			return size;
		};

		LOG(Level::TRACE, "Reading dataset 'I' from file '{}'", name);
		std::vector<double> buffer_i;
		const auto size = read_dataset("I", buffer_i);

		std::vector<double> buffer_q;
		LOG(Level::TRACE, "Reading dataset 'Q' from file '{}'", name);
		if (read_dataset("Q", buffer_q) != size)
		{
			LOG(Level::FATAL, "Dataset 'Q' is not the same size as dataset 'I' in file '{}'", name);
			throw std::runtime_error(R"(Dataset "Q" is not the same size as dataset "I" in file )" + name);
		}

		data.resize(size);
		for (size_t i = 0; i < size; ++i)
		{
			data[i] = ComplexType(buffer_i[i], buffer_q[i]);
		}
		LOG(Level::TRACE, "Read dataset successfully");
	}

	void addChunkToFile(HighFive::File& file, const std::vector<ComplexType>& data, const RealType time,
						const RealType fullscale, const unsigned count)
	{
		const unsigned size = data.size();

		const std::string base_chunk_name = "chunk_" + std::format("{:06}", count);
		const std::string i_chunk_name = base_chunk_name + "_I";
		const std::string q_chunk_name = base_chunk_name + "_Q";

		std::vector<RealType> i(size), q(size);
		std::ranges::transform(data, i.begin(), [](const ComplexType& c) { return c.real(); });
		std::ranges::transform(data, q.begin(), [](const ComplexType& c) { return c.imag(); });

		auto write_chunk = [&](const std::string& chunkName, const std::vector<RealType>& chunkData)
		{
			try
			{
				HighFive::DataSet dataset =
					file.createDataSet<RealType>(chunkName, HighFive::DataSpace::From(chunkData));
				dataset.write(chunkData);
			}
			catch (const HighFive::Exception& err)
			{
				LOG(Level::FATAL, "Error while writing data to HDF5 file: {}", err.what());
				throw std::runtime_error("Error while writing data to HDF5 file: " + chunkName + " - " + err.what());
			}
		};

		auto set_chunk_attributes = [&](const std::string& chunkName)
		{
			try
			{
				HighFive::DataSet dataset = file.getDataSet(chunkName);
				dataset.createAttribute("time", time);
				dataset.createAttribute("rate", params::rate());
				dataset.createAttribute("fullscale", fullscale);
			}
			catch (const HighFive::Exception& err)
			{
				LOG(Level::FATAL, "Error while setting attributes on chunk: {}", err.what());
				throw std::runtime_error("Error while setting attributes on chunk: " + chunkName + " - " + err.what());
			}
		};

		write_chunk(i_chunk_name, i);
		write_chunk(q_chunk_name, q);

		set_chunk_attributes(i_chunk_name);
		set_chunk_attributes(q_chunk_name);
	}

	std::vector<std::vector<RealType>> readPattern(const std::string& name, const std::string& datasetName)
	{
		try
		{
			LOG(Level::TRACE, "Reading dataset '{}' from file '{}'", datasetName, name);
			const HighFive::File file(name, HighFive::File::ReadOnly);

			const auto dataset = file.getDataSet(datasetName);

			const auto dataspace = dataset.getSpace();
			const auto dims = dataspace.getDimensions();

			if (dims.size() != 2)
			{
				LOG(Level::FATAL, "Invalid dataset dimensions for '{}' in file '{}'", datasetName, name);
				throw std::runtime_error(
					std::format(R"(Invalid dataset dimensions for "{}" in file "{}")", datasetName, name));
			}

			LOG(Level::TRACE, "Reading dataset with dimensions {}x{}", dims[0], dims[1]);

			std::vector data(dims[0], std::vector<RealType>(dims[1]));
			dataset.read(data);

			LOG(Level::TRACE, "Read dataset successfully");

			return data;
		}
		catch (const HighFive::Exception& err)
		{
			LOG(Level::FATAL, "Error handling HDF5 file: {}", err.what());
			throw std::runtime_error("Error handling HDF5 file: " + std::string(err.what()));
		}
	}
}
