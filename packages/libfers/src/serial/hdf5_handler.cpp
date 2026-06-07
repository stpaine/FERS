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
#include <string>
#include <vector>

#include "core/logging.h"
#include "core/parameters.h"

using logging::Level;

namespace serial
{
	std::mutex hdf5_global_mutex;

	namespace
	{
		template <typename T>
		void createOptionalAttribute(HighFive::File& file, const std::string& name, const std::optional<T>& value)
		{
			if (value.has_value())
			{
				file.createAttribute(name, *value);
			}
		}

		template <typename T>
		void createOptionalUnsignedAttribute(HighFive::File& file, const std::string& name,
											 const std::optional<T>& value)
		{
			if (value.has_value())
			{
				file.createAttribute(name, static_cast<unsigned long long>(*value));
			}
		}

		void writeBaseMetadataAttributes(HighFive::File& file, const core::OutputFileMetadata& metadata)
		{
			file.createAttribute("fers_metadata_schema_version", 1U);
			file.createAttribute("fers_metadata_json", core::outputFileMetadataToJsonString(metadata));
			file.createAttribute("receiver_id", static_cast<unsigned long long>(metadata.receiver_id));
			file.createAttribute("receiver_name", metadata.receiver_name);
			file.createAttribute("data_mode", metadata.mode);
			if (metadata.sampling_rate > 0.0)
			{
				file.createAttribute("output_sampling_rate", metadata.sampling_rate);
			}
			file.createAttribute("total_samples", static_cast<unsigned long long>(metadata.total_samples));
			file.createAttribute("sample_start", static_cast<unsigned long long>(metadata.sample_start));
			file.createAttribute("sample_end_exclusive",
								 static_cast<unsigned long long>(metadata.sample_end_exclusive));
			file.createAttribute("streaming_segment_count",
								 static_cast<unsigned long long>(metadata.streaming_segments.size()));
			file.createAttribute("fmcw_source_count", static_cast<unsigned long long>(metadata.fmcw_sources.size()));
			file.createAttribute("fmcw_dechirp_mode", metadata.fmcw_dechirp_mode);
			file.createAttribute("fmcw_dechirp_reference_source", metadata.fmcw_dechirp_reference_source);
			file.createAttribute("fmcw_if_decimation_enabled", metadata.fmcw_if_decimation_enabled);
			file.createAttribute("fmcw_if_legacy_full_rate", metadata.fmcw_if_legacy_full_rate);
			file.createAttribute("fmcw_if_group_delay_compensated", metadata.fmcw_if_group_delay_compensated);
		}

		void writeFmcwIfAttributes(HighFive::File& file, const core::OutputFileMetadata& metadata)
		{
			createOptionalAttribute(file, "fmcw_if_requested_sample_rate", metadata.fmcw_if_requested_sample_rate);
			createOptionalAttribute(file, "fmcw_if_sample_rate", metadata.fmcw_if_sample_rate);
			createOptionalAttribute(file, "fmcw_if_input_sample_rate", metadata.fmcw_if_input_sample_rate);
			createOptionalUnsignedAttribute(file, "fmcw_if_resample_numerator", metadata.fmcw_if_resample_numerator);
			createOptionalUnsignedAttribute(file, "fmcw_if_resample_denominator",
											metadata.fmcw_if_resample_denominator);
			createOptionalAttribute(file, "fmcw_if_decimation_factor", metadata.fmcw_if_decimation_factor);
			createOptionalAttribute(file, "fmcw_if_filter_bandwidth", metadata.fmcw_if_filter_bandwidth);
			createOptionalAttribute(file, "fmcw_if_filter_transition_width", metadata.fmcw_if_filter_transition_width);
			createOptionalAttribute(file, "fmcw_if_filter_stopband", metadata.fmcw_if_filter_stopband);
			createOptionalAttribute(file, "fmcw_if_filter_group_delay_seconds",
									metadata.fmcw_if_filter_group_delay_seconds);
			createOptionalUnsignedAttribute(file, "fmcw_if_compensated_integer_delay_samples",
											metadata.fmcw_if_compensated_integer_delay_samples);
			createOptionalAttribute(file, "fmcw_if_compensated_fractional_delay_samples",
									metadata.fmcw_if_compensated_fractional_delay_samples);
			createOptionalUnsignedAttribute(file, "fmcw_if_warmup_discard_samples",
											metadata.fmcw_if_warmup_discard_samples);
			createOptionalUnsignedAttribute(file, "fmcw_if_phase_refinement", metadata.fmcw_if_phase_refinement);
			createOptionalAttribute(file, "fmcw_if_timing_error_seconds", metadata.fmcw_if_timing_error_seconds);
			createOptionalAttribute(file, "fmcw_if_phase_error_radians", metadata.fmcw_if_phase_error_radians);
			createOptionalAttribute(file, "fmcw_if_noise_variance", metadata.fmcw_if_noise_variance);
		}

		void writeFmcwWaveformAttributes(HighFive::File& file, const std::string& prefix,
										 const core::FmcwMetadata& waveform, const bool write_counts)
		{
			file.createAttribute(prefix + "waveform_shape", waveform.waveform_shape);
			file.createAttribute(prefix + "chirp_bandwidth", waveform.chirp_bandwidth);
			file.createAttribute(prefix + "chirp_duration", waveform.chirp_duration);
			file.createAttribute(prefix + "chirp_rate", waveform.chirp_rate);
			file.createAttribute(prefix + "start_frequency_offset", waveform.start_frequency_offset);
			if (waveform.waveform_shape == "linear")
			{
				file.createAttribute(prefix + "chirp_period", waveform.chirp_period);
				file.createAttribute(prefix + "chirp_direction", waveform.chirp_direction);
				if (write_counts)
				{
					file.createAttribute(prefix + "chirp_rate_signed", waveform.chirp_rate_signed);
					createOptionalUnsignedAttribute(file, prefix + "chirp_count", waveform.chirp_count);
				}
			}
			else if (waveform.waveform_shape == "triangle" && waveform.triangle_period.has_value())
			{
				file.createAttribute(prefix + "triangle_period", *waveform.triangle_period);
				if (write_counts)
				{
					createOptionalUnsignedAttribute(file, prefix + "triangle_count", waveform.triangle_count);
				}
			}
		}

		void writeDechirpReferenceAttributes(HighFive::File& file, const core::OutputFileMetadata& metadata)
		{
			createOptionalUnsignedAttribute(file, "fmcw_dechirp_reference_transmitter_id",
											metadata.fmcw_dechirp_reference_transmitter_id);
			createOptionalAttribute(file, "fmcw_dechirp_reference_transmitter_name",
									metadata.fmcw_dechirp_reference_transmitter_name);
			createOptionalUnsignedAttribute(file, "fmcw_dechirp_reference_waveform_id",
											metadata.fmcw_dechirp_reference_waveform_id);
			createOptionalAttribute(file, "fmcw_dechirp_reference_waveform_name",
									metadata.fmcw_dechirp_reference_waveform_name);
			if (metadata.fmcw_dechirp_reference_waveform.has_value())
			{
				writeFmcwWaveformAttributes(file, "fmcw_dechirp_reference_", *metadata.fmcw_dechirp_reference_waveform,
											false);
			}
		}

		void writeStreamingSegmentFmcwAttributes(HighFive::File& file, const core::OutputFileMetadata& metadata)
		{
			std::vector<RealType> streaming_first_chirp_starts;
			std::vector<unsigned long long> streaming_emitted_chirp_counts;
			std::vector<RealType> streaming_first_triangle_starts;
			std::vector<unsigned long long> streaming_emitted_triangle_counts;
			for (const auto& segment : metadata.streaming_segments)
			{
				if (segment.first_chirp_start_time.has_value())
				{
					streaming_first_chirp_starts.push_back(*segment.first_chirp_start_time);
				}
				if (segment.emitted_chirp_count.has_value())
				{
					streaming_emitted_chirp_counts.push_back(
						static_cast<unsigned long long>(*segment.emitted_chirp_count));
				}
				if (segment.first_triangle_start_time.has_value())
				{
					streaming_first_triangle_starts.push_back(*segment.first_triangle_start_time);
				}
				if (segment.emitted_triangle_count.has_value())
				{
					streaming_emitted_triangle_counts.push_back(
						static_cast<unsigned long long>(*segment.emitted_triangle_count));
				}
			}
			if (!streaming_first_chirp_starts.empty())
			{
				auto attr = file.createAttribute<RealType>("streaming_first_chirp_start_time",
														   HighFive::DataSpace::From(streaming_first_chirp_starts));
				attr.write(streaming_first_chirp_starts);
			}
			if (!streaming_emitted_chirp_counts.empty())
			{
				auto attr = file.createAttribute<unsigned long long>(
					"streaming_emitted_chirp_count", HighFive::DataSpace::From(streaming_emitted_chirp_counts));
				attr.write(streaming_emitted_chirp_counts);
			}
			if (!streaming_first_triangle_starts.empty())
			{
				auto attr = file.createAttribute<RealType>("streaming_first_triangle_start_time",
														   HighFive::DataSpace::From(streaming_first_triangle_starts));
				attr.write(streaming_first_triangle_starts);
			}
			if (!streaming_emitted_triangle_counts.empty())
			{
				auto attr = file.createAttribute<unsigned long long>(
					"streaming_emitted_triangle_count", HighFive::DataSpace::From(streaming_emitted_triangle_counts));
				attr.write(streaming_emitted_triangle_counts);
			}
		}

		void writeFmcwAttributes(HighFive::File& file, const core::OutputFileMetadata& metadata)
		{
			if (!metadata.fmcw.has_value())
			{
				return;
			}
			writeFmcwWaveformAttributes(file, "fmcw_", *metadata.fmcw, true);
			writeStreamingSegmentFmcwAttributes(file, metadata);
		}
	}

	void writeOutputFileMetadataAttributes(HighFive::File& file, const core::OutputFileMetadata& metadata)
	{
		writeBaseMetadataAttributes(file, metadata);
		writeFmcwIfAttributes(file, metadata);
		writeDechirpReferenceAttributes(file, metadata);
		writeFmcwAttributes(file, metadata);
	}

	void readPulseData(const std::string& name, std::vector<ComplexType>& data)
	{
		std::scoped_lock const lock(hdf5_global_mutex);

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
						const RealType fullscale, const unsigned count, const core::PulseChunkMetadata* metadata)
	{
		std::scoped_lock const lock(hdf5_global_mutex);

		const std::size_t size = data.size();

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
				if (metadata != nullptr)
				{
					dataset.createAttribute("chunk_index", metadata->chunk_index);
					dataset.createAttribute("sample_count", static_cast<unsigned long long>(metadata->sample_count));
					dataset.createAttribute("sample_start", static_cast<unsigned long long>(metadata->sample_start));
					dataset.createAttribute("sample_end_exclusive",
											static_cast<unsigned long long>(metadata->sample_end_exclusive));
				}
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
		std::scoped_lock const lock(hdf5_global_mutex);
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
