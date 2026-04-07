// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "finalizer_pipeline.h"

#include <algorithm>
#include <cmath>
#include <highfive/highfive.hpp>
#include <ranges>

#include "core/logging.h"
#include "core/parameters.h"
#include "processing/signal_processor.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/hdf5_handler.h"
#include "serial/response.h"
#include "signal/dsp_filters.h"
#include "simulation/channel_model.h"
#include "timing/timing.h"

namespace processing::pipeline
{
	void advanceTimingModel(timing::Timing* timing_model, const radar::Receiver* receiver, const RealType rate)
	{
		if ((timing_model == nullptr) || !timing_model->isEnabled())
		{
			return;
		}

		if (timing_model->getSyncOnPulse())
		{
			timing_model->reset();
			timing_model->skipSamples(static_cast<long>(std::floor(rate * receiver->getWindowSkip())));
		}
		else
		{
			const RealType inter_pulse_skip_duration = 1.0 / receiver->getWindowPrf() - receiver->getWindowLength();
			const auto samples_to_skip = static_cast<long>(std::floor(rate * inter_pulse_skip_duration));
			timing_model->skipSamples(samples_to_skip);
		}
	}

	std::tuple<RealType, RealType> calculateJitteredStart(const RealType ideal_start, const RealType first_phase_noise,
														  const RealType carrier_freq, const RealType rate)
	{
		const RealType actual_start = ideal_start + first_phase_noise / (2.0 * PI * carrier_freq);
		const RealType rounded_start = std::round(actual_start * rate) / rate;
		const RealType fractional_delay = actual_start * rate - std::round(actual_start * rate);
		return {rounded_start, fractional_delay};
	}

	void applyCwInterference(std::span<ComplexType> window, const RealType actual_start, const RealType dt,
							 const radar::Receiver* receiver, const std::vector<radar::Transmitter*>& cw_sources,
							 const std::vector<std::unique_ptr<radar::Target>>* targets)
	{
		RealType t_sample = actual_start;
		for (auto& window_sample : window)
		{
			ComplexType cw_interference_sample{0.0, 0.0};
			for (const auto* cw_source : cw_sources)
			{
				if (!receiver->checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT))
				{
					cw_interference_sample +=
						simulation::calculateDirectPathContribution(cw_source, receiver, t_sample);
				}
				for (const auto& target_ptr : *targets)
				{
					cw_interference_sample +=
						simulation::calculateReflectedPathContribution(cw_source, receiver, target_ptr.get(), t_sample);
				}
			}
			window_sample += cw_interference_sample;
			t_sample += dt;
		}
	}

	void applyPulsedInterference(std::vector<ComplexType>& iq_buffer,
								 const std::vector<std::unique_ptr<serial::Response>>& interference_log)
	{
		for (const auto& response : interference_log)
		{
			unsigned psize;
			RealType prate;
			const auto rendered_pulse = response->renderBinary(prate, psize, 0.0);

			const RealType dt_sim = 1.0 / prate;
			const auto start_index = static_cast<size_t>((response->startTime() - params::startTime()) / dt_sim);

			for (size_t i = 0; i < psize; ++i)
			{
				if (start_index + i < iq_buffer.size())
				{
					iq_buffer[start_index + i] += rendered_pulse[i];
				}
			}
		}
	}

	void addPhaseNoiseToWindow(std::span<const RealType> noise, std::span<ComplexType> window)
	{
		for (auto [n, w] : std::views::zip(noise, window))
		{
			w *= std::polar(1.0, n);
		}
	}

	RealType applyDownsamplingAndQuantization(std::vector<ComplexType>& buffer)
	{
		if (params::oversampleRatio() > 1)
		{
			buffer = std::move(fers_signal::downsample(buffer));
		}
		return quantizeAndScaleWindow(buffer);
	}

	void exportCwToHdf5(const std::string& filename, const std::vector<ComplexType>& iq_buffer,
						const RealType fullscale, const RealType ref_freq)
	{
		std::scoped_lock lock(serial::hdf5_global_mutex);
		try
		{
			HighFive::File file(filename, HighFive::File::Truncate);

			std::vector<RealType> i_data(iq_buffer.size());
			std::vector<RealType> q_data(iq_buffer.size());
			std::ranges::transform(iq_buffer, i_data.begin(), [](const auto& c) { return c.real(); });
			std::ranges::transform(iq_buffer, q_data.begin(), [](const auto& c) { return c.imag(); });

			HighFive::DataSet i_dataset = file.createDataSet<RealType>("I_data", HighFive::DataSpace::From(i_data));
			i_dataset.write(i_data);
			HighFive::DataSet q_dataset = file.createDataSet<RealType>("Q_data", HighFive::DataSpace::From(q_data));
			q_dataset.write(q_data);

			file.createAttribute("sampling_rate", params::rate());
			file.createAttribute("start_time", params::startTime());
			file.createAttribute("fullscale", fullscale);
			file.createAttribute("reference_carrier_frequency", ref_freq);

			LOG(logging::Level::INFO, "Successfully exported CW data to '{}'", filename);
		}
		catch (const HighFive::Exception& err)
		{
			LOG(logging::Level::FATAL, "Error writing CW data to HDF5 file '{}': {}", filename, err.what());
		}
	}
}
