// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

#include <algorithm>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <nlohmann/json.hpp>
#include <random>
#include <string>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/parameters.h"
#include "core/sim_id.h"
#include "core/world.h"
#include "math/path.h"
#include "math/rotation_path.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/json_serializer.h"
#include "signal/radar_signal.h"
#include "timing/prototype_timing.h"
#include "timing/timing.h"

using json = nlohmann::json;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		ParamGuard(const ParamGuard&) = delete;
		ParamGuard& operator=(const ParamGuard&) = delete;
		ParamGuard(ParamGuard&&) = delete;
		ParamGuard& operator=(ParamGuard&&) = delete;
		~ParamGuard() { params::params = saved; }
	};

	// Helper to recursively sort JSON arrays by "id" to ensure deterministic
	// comparison, as core::World uses unordered_maps internally.
	void sortJsonArrays(json& j)
	{
		if (j.is_object())
		{
			for (const auto& [key, val] : j.items())
			{
				sortJsonArrays(val);
			}
		}
		else if (j.is_array())
		{
			if (!j.empty() && j[0].is_object() && j[0].contains("id"))
			{
				// NOLINTNEXTLINE(modernize-use-ranges): nlohmann::json iterators fail ranges::sort concepts here.
				std::sort(j.begin(), j.end(),
						  [](const json& a, const json& b)
						  {
							  // Parse stringified uint64_t IDs safely for accurate numerical sorting
							  if (a["id"].is_string() && b["id"].is_string())
								  return std::stoull(a["id"].get<std::string>()) <
									  std::stoull(b["id"].get<std::string>());
							  if (a["id"].is_number() && b["id"].is_number())
								  return a["id"].get<uint64_t>() < b["id"].get<uint64_t>();
							  return false;
						  });
			}
			for (auto& val : j)
			{
				sortJsonArrays(val);
			}
		}
	}

	// Smart recursive JSON comparator that handles floating point drift
	// and prints the exact path of the failure.
	bool compareJson(const json& a, const json& b, const std::string& path = "");

	void reportMissingObjectKeys(const json& a, const json& b, const std::string& path)
	{
		for (const auto& [k, v] : a.items())
		{
			if (!b.contains(k))
			{
				UNSCOPED_INFO("Key missing in B: " << path << "/" << k);
			}
		}
		for (const auto& [k, v] : b.items())
		{
			if (!a.contains(k))
			{
				UNSCOPED_INFO("Key missing in A: " << path << "/" << k);
			}
		}
	}

	bool compareJsonObject(const json& a, const json& b, const std::string& path)
	{
		if (a.size() != b.size())
		{
			UNSCOPED_INFO("Object size mismatch at " << path << ": " << a.size() << " vs " << b.size());
			reportMissingObjectKeys(a, b, path);
			return false;
		}
		for (const auto& [key, val] : a.items())
		{
			std::string child_path = path;
			child_path += '/';
			child_path += key;
			if (!compareJson(val, b[key], child_path))
			{
				return false;
			}
		}
		return true;
	}

	bool compareJsonArray(const json& a, const json& b, const std::string& path)
	{
		if (a.size() != b.size())
		{
			UNSCOPED_INFO("Array size mismatch at " << path << ": " << a.size() << " vs " << b.size());
			return false;
		}
		for (size_t i = 0; i < a.size(); ++i)
		{
			if (!compareJson(a[i], b[i], path + "[" + std::to_string(i) + "]"))
			{
				return false;
			}
		}
		return true;
	}

	bool compareJsonFloat(const json& a, const json& b, const std::string& path)
	{
		const double va = a.get<double>();
		const double vb = b.get<double>();
		// Allow a small epsilon for float serialization round-tripping
		if (std::abs(va - vb) > 1e-5 && std::abs(va - vb) / std::max(std::abs(va), std::abs(vb)) > 1e-5)
		{
			UNSCOPED_INFO("Float mismatch at " << path << ": " << va << " vs " << vb);
			return false;
		}
		return true;
	}

	bool compareJson(const json& a, const json& b, const std::string& path)
	{
		if (a.type() != b.type())
		{
			UNSCOPED_INFO("Type mismatch at " << path << ": " << a.type_name() << " vs " << b.type_name());
			return false;
		}
		if (a.is_object())
		{
			return compareJsonObject(a, b, path);
		}
		if (a.is_array())
		{
			return compareJsonArray(a, b, path);
		}
		if (a.is_number_float())
		{
			return compareJsonFloat(a, b, path);
		}

		if (a != b)
		{
			UNSCOPED_INFO("Value mismatch at " << path << ": " << a << " vs " << b);
			return false;
		}
		return true;
	}

	void addStressAssets(core::World& world, std::mt19937& rng, std::vector<SimId>& wave_ids,
						 std::vector<SimId>& ant_ids, std::vector<SimId>& time_ids)
	{
		std::uniform_real_distribution<RealType> dist_real(0.1, 1000.0);

		for (size_t i = 0; i < 20; ++i)
		{
			SimId const w_id = SimIdGenerator::instance().generateId(ObjectType::Waveform);
			auto sig = std::make_unique<fers_signal::CwSignal>();
			auto wave = std::make_unique<fers_signal::RadarSignal>("wave_" + std::to_string(i), dist_real(rng),
																   1e9 + dist_real(rng), 1.0, std::move(sig), w_id);
			world.add(std::move(wave));
			wave_ids.push_back(w_id);

			SimId const a_id = SimIdGenerator::instance().generateId(ObjectType::Antenna);
			std::unique_ptr<antenna::Antenna> ant;
			switch (i % 5)
			{
			case 0:
				ant = std::make_unique<antenna::Isotropic>("ant_" + std::to_string(i), a_id);
				break;
			case 1:
				ant = std::make_unique<antenna::Sinc>("ant_" + std::to_string(i), 1.0, 2.0, 3.0, a_id);
				break;
			case 2:
				ant = std::make_unique<antenna::Gaussian>("ant_" + std::to_string(i), 1.5, 1.5, a_id);
				break;
			case 3:
				ant = std::make_unique<antenna::SquareHorn>("ant_" + std::to_string(i), 0.5, a_id);
				break;
			case 4:
				ant = std::make_unique<antenna::Parabolic>("ant_" + std::to_string(i), 2.0, a_id);
				break;
			}
			ant->setEfficiencyFactor(0.95);
			world.add(std::move(ant));
			ant_ids.push_back(a_id);

			SimId const t_id = SimIdGenerator::instance().generateId(ObjectType::Timing);
			auto tim = std::make_unique<timing::PrototypeTiming>("time_" + std::to_string(i), t_id);
			tim->setFrequency(10e6);
			tim->setFreqOffset(dist_real(rng));
			tim->setRandomFreqOffsetStdev(dist_real(rng));
			tim->setPhaseOffset(dist_real(rng));
			tim->setRandomPhaseOffsetStdev(dist_real(rng));
			tim->setAlpha(1.0, 0.5);
			tim->setAlpha(2.0, 0.25);
			if (i % 2 == 0)
			{
				tim->setSyncOnPulse();
			}
			world.add(std::move(tim));
			time_ids.push_back(t_id);
		}
	}

	void addStressPlatform(core::World& world, const size_t i, std::mt19937& rng,
						   std::uniform_real_distribution<RealType>& dist_real,
						   std::uniform_int_distribution<unsigned>& seed_dist, const std::vector<SimId>& wave_ids,
						   const std::vector<SimId>& ant_ids, const std::vector<SimId>& time_ids)
	{
		SimId const p_id = SimIdGenerator::instance().generateId(ObjectType::Platform);
		auto plat = std::make_unique<radar::Platform>("platform_" + std::to_string(i), p_id);

		plat->getMotionPath()->setInterp(math::Path::InterpType::INTERP_CUBIC);
		for (size_t wp = 0; wp < 10; ++wp)
		{
			plat->getMotionPath()->addCoord(
				{math::Vec3(dist_real(rng), dist_real(rng), dist_real(rng)), static_cast<RealType>(wp) * 0.01});
		}
		plat->getMotionPath()->finalize();

		if (i % 2 == 0)
		{
			plat->getRotationPath()->setInterp(math::RotationPath::InterpType::INTERP_LINEAR);
			for (size_t wp = 0; wp < 10; ++wp)
			{
				plat->getRotationPath()->addCoord(
					{dist_real(rng) * (PI / 180.0), dist_real(rng) * (PI / 180.0), static_cast<RealType>(wp) * 0.01});
			}
		}
		else
		{
			plat->getRotationPath()->setInterp(math::RotationPath::InterpType::INTERP_CONSTANT);
			math::RotationCoord const start{dist_real(rng) * (PI / 180.0), dist_real(rng) * (PI / 180.0), 0.0};
			math::RotationCoord const rate{dist_real(rng) * (PI / 180.0), dist_real(rng) * (PI / 180.0), 0.0};
			plat->getRotationPath()->setConstantRate(start, rate);
		}
		plat->getRotationPath()->finalize();

		auto* proto_tim = world.findTiming(time_ids[i % time_ids.size()]);

		auto tx =
			std::make_unique<radar::Transmitter>(plat.get(), "tx_" + std::to_string(i), radar::OperationMode::CW_MODE);
		tx->setWave(world.findWaveform(wave_ids[i % wave_ids.size()]));
		tx->setAntenna(world.findAntenna(ant_ids[i % ant_ids.size()]));
		auto tx_tim = std::make_shared<timing::Timing>(proto_tim->getName(), seed_dist(rng), proto_tim->getId());
		tx_tim->initializeModel(proto_tim);
		tx->setTiming(tx_tim);
		tx->setSchedule({{0.01, 0.04}, {0.06, 0.09}});
		world.add(std::move(tx));

		auto rx = std::make_unique<radar::Receiver>(plat.get(), "rx_" + std::to_string(i), seed_dist(rng),
													radar::OperationMode::CW_MODE);
		rx->setAntenna(world.findAntenna(ant_ids[(i + 1) % ant_ids.size()]));
		auto rx_tim = std::make_shared<timing::Timing>(proto_tim->getName(), seed_dist(rng), proto_tim->getId());
		rx_tim->initializeModel(proto_tim);
		rx->setTiming(rx_tim);
		rx->setNoiseTemperature(290.0);
		if (i % 2 == 0)
		{
			rx->setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
		}
		if (i % 3 == 0)
		{
			rx->setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
		}
		rx->setSchedule({{0.02, 0.08}});
		world.add(std::move(rx));

		auto mono_tx = std::make_unique<radar::Transmitter>(plat.get(), "mono_tx_" + std::to_string(i),
															radar::OperationMode::CW_MODE);
		mono_tx->setWave(world.findWaveform(wave_ids[(i + 2) % wave_ids.size()]));
		mono_tx->setAntenna(world.findAntenna(ant_ids[(i + 2) % ant_ids.size()]));
		auto mono_tx_tim = std::make_shared<timing::Timing>(proto_tim->getName(), seed_dist(rng), proto_tim->getId());
		mono_tx_tim->initializeModel(proto_tim);
		mono_tx->setTiming(mono_tx_tim);
		mono_tx->setSchedule({{0.01, 0.09}});

		auto mono_rx = std::make_unique<radar::Receiver>(plat.get(), "mono_rx_" + std::to_string(i), seed_dist(rng),
														 radar::OperationMode::CW_MODE);
		mono_rx->setAntenna(world.findAntenna(ant_ids[(i + 2) % ant_ids.size()]));
		auto mono_rx_tim = std::make_shared<timing::Timing>(proto_tim->getName(), seed_dist(rng), proto_tim->getId());
		mono_rx_tim->initializeModel(proto_tim);
		mono_rx->setTiming(mono_rx_tim);
		mono_rx->setNoiseTemperature(300.0);
		mono_rx->setSchedule({{0.01, 0.09}});

		mono_tx->setAttached(mono_rx.get());
		mono_rx->setAttached(mono_tx.get());
		world.add(std::move(mono_tx));
		world.add(std::move(mono_rx));

		auto tgt = radar::createIsoTarget(plat.get(), "tgt_" + std::to_string(i), dist_real(rng), seed_dist(rng));
		if (i % 2 == 0)
		{
			tgt->setFluctuationModel(std::make_unique<radar::RcsChiSquare>(tgt->getRngEngine(), 2.0));
		}
		else
		{
			tgt->setFluctuationModel(std::make_unique<radar::RcsConst>());
		}
		world.add(std::move(tgt));
		world.add(std::move(plat));
	}

	void buildStressWorld(core::World& world, size_t num_platforms, std::mt19937& rng)
	{
		std::uniform_real_distribution<RealType> dist_real(0.1, 1000.0);
		std::uniform_int_distribution<unsigned> seed_dist;
		std::vector<SimId> wave_ids, ant_ids, time_ids;

		addStressAssets(world, rng, wave_ids, ant_ids, time_ids);
		for (size_t i = 0; i < num_platforms; ++i)
		{
			addStressPlatform(world, i, rng, dist_real, seed_dist, wave_ids, ant_ids, time_ids);
		}
	}
}

TEST_CASE("JSON Serializer Stress Test and Round-Trip Validation", "[serial][json][stress]")
{
	ParamGuard const guard;
	params::params.reset();

	// Keep time and rate low to prevent massive CW buffer allocations during json_to_world
	params::setTime(0.0, 0.1);
	params::setRate(1000.0);
	params::setCoordinateSystem(params::CoordinateFrame::UTM, 34, false);
	params::params.simulation_name = "StressTestSim";
	params::params.random_seed = 1337;

	std::mt19937 master_seeder(1337);

	core::World original_world;

	// 10 Platforms is plenty to stress the JSON parser without taking forever.
	// This will generate 10 standalone Txs, 10 standalone Rxs, 10 Monostatic pairs (20 objects), and 10 Targets.
	constexpr size_t STRESS_PLATFORM_COUNT = 10;

	buildStressWorld(original_world, STRESS_PLATFORM_COUNT, master_seeder);

	json j_original;
	json j_roundtrip;

	SECTION("Deep Equality Round-Trip Validation")
	{
		j_original = serial::world_to_json(original_world);
		REQUIRE(!j_original.empty());

		core::World new_world;
		std::mt19937 new_seeder(1337);

		// Deserialize
		REQUIRE_NOTHROW(serial::json_to_world(j_original, new_world, new_seeder));

		// Verify object counts
		REQUIRE(new_world.getPlatforms().size() == STRESS_PLATFORM_COUNT);
		REQUIRE(new_world.getTransmitters().size() == STRESS_PLATFORM_COUNT * 2); // Standalone + Monostatic Tx
		REQUIRE(new_world.getReceivers().size() == STRESS_PLATFORM_COUNT * 2); // Standalone + Monostatic Rx
		REQUIRE(new_world.getTargets().size() == STRESS_PLATFORM_COUNT);
		REQUIRE(new_world.getWaveforms().size() == 20);
		REQUIRE(new_world.getAntennas().size() == 20);
		REQUIRE(new_world.getTimings().size() == 20);

		// Serialize again
		j_roundtrip = serial::world_to_json(new_world);
		REQUIRE(!j_roundtrip.empty());

		// Sort arrays by ID to ensure deterministic comparison
		sortJsonArrays(j_original);
		sortJsonArrays(j_roundtrip);

		// Use the smart comparator instead of == to handle float drift
		// and to print exact JSON paths on failure.
		bool const is_match = compareJson(j_original, j_roundtrip, "/root");
		REQUIRE(is_match);
	}

	SECTION("Performance Benchmarks")
	{
		BENCHMARK("world_to_json (Serialization)") { return serial::world_to_json(original_world); };

		json j_bench = serial::world_to_json(original_world);

		BENCHMARK("json_to_world (Deserialization)")
		{
			core::World bench_world;
			std::mt19937 bench_seeder(42);
			serial::json_to_world(j_bench, bench_world, bench_seeder);
			return bench_world.getPlatforms().size();
		};
	}
}
