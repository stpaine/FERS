// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "finalizer.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <format>
#include <highfive/highfive.hpp>
#include <vector>

#include "core/logging.h"
#include "core/parameters.h"
#include "core/rendering_job.h"
#include "core/sim_threading.h"
#include "processing/finalizer_pipeline.h"
#include "processing/signal_processor.h"
#include "radar/receiver.h"
#include "serial/hdf5_handler.h"
#include "timing/timing.h"

namespace processing
{
	void runPulsedFinalizer(radar::Receiver* receiver, const std::vector<std::unique_ptr<radar::Target>>* targets,
							std::shared_ptr<core::ProgressReporter> reporter)
	{
		const auto timing_model = receiver->getTiming()->clone();
		if (!timing_model)
		{
			LOG(logging::Level::FATAL, "Failed to clone timing model for receiver '{}'", receiver->getName());
			return;
		}

		const auto hdf5_filename = std::format("{}_results.h5", receiver->getName());

		std::unique_ptr<HighFive::File> h5_file;
		{
			std::scoped_lock lock(serial::hdf5_global_mutex);
			h5_file = std::make_unique<HighFive::File>(hdf5_filename, HighFive::File::Truncate);
		}

		unsigned chunk_index = 0;

		LOG(logging::Level::INFO, "Finalizer thread started for receiver '{}'. Outputting to '{}'.",
			receiver->getName(), hdf5_filename);

		auto last_report_time = std::chrono::steady_clock::now();
		const auto report_interval = std::chrono::milliseconds(100);
		const RealType rate = params::rate() * params::oversampleRatio();
		const RealType dt = 1.0 / rate;

		while (true)
		{
			core::RenderingJob job;
			if (!receiver->waitAndDequeueFinalizerJob(job))
			{
				break; // Shutdown signal received
			}

			const auto window_samples = static_cast<unsigned>(std::ceil(job.duration * rate));
			std::vector pnoise(window_samples, 0.0);

			RealType actual_start = job.ideal_start_time;
			RealType frac_delay = 0.0;

			if (timing_model->isEnabled())
			{
				pipeline::advanceTimingModel(timing_model.get(), receiver, rate);
				std::ranges::generate(pnoise, [&] { return timing_model->getNextSample(); });
				std::tie(actual_start, frac_delay) = pipeline::calculateJitteredStart(
					job.ideal_start_time, pnoise[0], timing_model->getFrequency(), rate);
			}

			std::vector<ComplexType> window_buffer(window_samples);

			applyThermalNoise(window_buffer, receiver->getNoiseTemperature(receiver->getRotation(actual_start)),
							  receiver->getRngEngine());

			pipeline::applyCwInterference(window_buffer, actual_start, dt, receiver, job.active_cw_sources, targets);

			renderWindow(window_buffer, job.duration, actual_start, frac_delay, job.responses);

			if (timing_model->isEnabled())
			{
				pipeline::addPhaseNoiseToWindow(pnoise, window_buffer);
			}

			const RealType fullscale = pipeline::applyDownsamplingAndQuantization(window_buffer);

			serial::addChunkToFile(*h5_file, window_buffer, actual_start, fullscale, chunk_index++);

			if (reporter)
			{
				const auto now = std::chrono::steady_clock::now();
				if ((now - last_report_time) >= report_interval)
				{
					reporter->report(std::format("Exporting {}: Chunk {}", receiver->getName(), chunk_index),
									 static_cast<int>(chunk_index), 0);
					last_report_time = now;
				}
			}
		}

		{
			// Safe destruction of the HDF5 object inside a lock
			std::scoped_lock lock(serial::hdf5_global_mutex);
			h5_file.reset();
		}

		if (reporter)
		{
			reporter->report(std::format("Finished Exporting {}", receiver->getName()), 100, 100);
		}
		LOG(logging::Level::INFO, "Finalizer thread for receiver '{}' finished.", receiver->getName());
	}

	void finalizeCwReceiver(radar::Receiver* receiver, pool::ThreadPool* /*pool*/,
							std::shared_ptr<core::ProgressReporter> reporter)
	{
		LOG(logging::Level::INFO, "Finalization task started for CW receiver '{}'.", receiver->getName());
		if (reporter)
		{
			reporter->report(std::format("Finalizing CW Receiver {}", receiver->getName()), 0, 100);
		}

		auto& iq_buffer = receiver->getMutableCwData();
		if (iq_buffer.empty())
		{
			LOG(logging::Level::INFO, "No CW data to finalize for receiver '{}'.", receiver->getName());
			return;
		}

		if (reporter)
		{
			reporter->report(std::format("Rendering Interference for {}", receiver->getName()), 25, 100);
		}
		pipeline::applyPulsedInterference(iq_buffer, receiver->getPulsedInterferenceLog());

		const auto timing_model = receiver->getTiming()->clone();
		if (!timing_model)
		{
			LOG(logging::Level::FATAL, "Failed to clone timing model for CW receiver '{}'", receiver->getName());
			return;
		}

		if (reporter)
		{
			reporter->report(std::format("Applying Noise for {}", receiver->getName()), 50, 100);
		}
		applyThermalNoise(iq_buffer, receiver->getNoiseTemperature(), receiver->getRngEngine());

		if (timing_model->isEnabled())
		{
			std::vector pnoise(iq_buffer.size(), 0.0);
			std::ranges::generate(pnoise, [&] { return timing_model->getNextSample(); });
			pipeline::addPhaseNoiseToWindow(pnoise, iq_buffer);
		}

		const RealType fullscale = pipeline::applyDownsamplingAndQuantization(iq_buffer);

		if (reporter)
		{
			reporter->report(std::format("Writing HDF5 for {}", receiver->getName()), 75, 100);
		}

		const auto hdf5_filename = std::format("{}_results.h5", receiver->getName());
		pipeline::exportCwToHdf5(hdf5_filename, iq_buffer, fullscale, timing_model->getFrequency());

		if (reporter)
		{
			reporter->report(std::format("Finalized {}", receiver->getName()), 100, 100);
		}
	}
}
