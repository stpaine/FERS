#include <catch2/catch_test_macros.hpp>
#include <cstdint>

#include "core/sim_id.h"

TEST_CASE("SimIdGenerator encodes type and counter", "[core][sim_id]")
{
	auto& generator = SimIdGenerator::instance();

	const SimId id1 = generator.generateId(ObjectType::Platform);
	const SimId id2 = generator.generateId(ObjectType::Platform);
	const SimId id3 = generator.generateId(ObjectType::Receiver);

	REQUIRE(SimIdGenerator::getType(id1) == ObjectType::Platform);
	REQUIRE(SimIdGenerator::getType(id2) == ObjectType::Platform);
	REQUIRE(SimIdGenerator::getType(id3) == ObjectType::Receiver);

	const uint64_t counter1 = SimIdGenerator::getCounter(id1);
	const uint64_t counter2 = SimIdGenerator::getCounter(id2);
	const uint64_t counter3 = SimIdGenerator::getCounter(id3);

	REQUIRE(counter2 == counter1 + 1);
	REQUIRE(counter3 == counter2 + 1);
}

TEST_CASE("SimIdGenerator produces debug IDs", "[core][sim_id]")
{
	auto& generator = SimIdGenerator::instance();

	const SimId debug_id = generator.generateDebugId();
	REQUIRE(SimIdGenerator::getType(debug_id) == ObjectType::Debug);
	REQUIRE(SimIdGenerator::getCounter(debug_id) > 0);
}

TEST_CASE("SimId helpers decode manually constructed ids", "[core][sim_id]")
{
	constexpr uint64_t counter_mask = 0x0000FFFFFFFFFFFF;
	const uint64_t sample_counter = 0x0000000000ABCDEF & counter_mask;

	const SimId platform_id = (static_cast<uint64_t>(ObjectType::Platform) << 48) | sample_counter;
	const SimId debug_id = (static_cast<uint64_t>(ObjectType::Debug) << 48) | sample_counter;

	REQUIRE(SimIdGenerator::getType(platform_id) == ObjectType::Platform);
	REQUIRE(SimIdGenerator::getCounter(platform_id) == sample_counter);
	REQUIRE(SimIdGenerator::getType(debug_id) == ObjectType::Debug);
	REQUIRE(SimIdGenerator::getCounter(debug_id) == sample_counter);
}

TEST_CASE("SimIdGenerator encodes all object types", "[core][sim_id]")
{
	auto& generator = SimIdGenerator::instance();

	const SimId unknown_id = generator.generateId(ObjectType::Unknown);
	const SimId platform_id = generator.generateId(ObjectType::Platform);
	const SimId transmitter_id = generator.generateId(ObjectType::Transmitter);
	const SimId receiver_id = generator.generateId(ObjectType::Receiver);
	const SimId target_id = generator.generateId(ObjectType::Target);
	const SimId antenna_id = generator.generateId(ObjectType::Antenna);
	const SimId waveform_id = generator.generateId(ObjectType::Waveform);
	const SimId timing_id = generator.generateId(ObjectType::Timing);

	REQUIRE(SimIdGenerator::getType(unknown_id) == ObjectType::Unknown);
	REQUIRE(SimIdGenerator::getType(platform_id) == ObjectType::Platform);
	REQUIRE(SimIdGenerator::getType(transmitter_id) == ObjectType::Transmitter);
	REQUIRE(SimIdGenerator::getType(receiver_id) == ObjectType::Receiver);
	REQUIRE(SimIdGenerator::getType(target_id) == ObjectType::Target);
	REQUIRE(SimIdGenerator::getType(antenna_id) == ObjectType::Antenna);
	REQUIRE(SimIdGenerator::getType(waveform_id) == ObjectType::Waveform);
	REQUIRE(SimIdGenerator::getType(timing_id) == ObjectType::Timing);
}

TEST_CASE("SimIdGenerator preserves counter mask", "[core][sim_id]")
{
	auto& generator = SimIdGenerator::instance();

	const SimId id = generator.generateId(ObjectType::Target);
	const uint64_t counter = SimIdGenerator::getCounter(id);
	constexpr uint64_t counter_mask = 0x0000FFFFFFFFFFFF;

	REQUIRE((counter & counter_mask) == counter);
}

TEST_CASE("SimIdGenerator edge cases not directly testable", "[core][sim_id]")
{
	// TODO: Cannot reliably trigger the 48-bit counter overflow in unit tests.
	// TODO: Cannot validate the debug-type assert without build-time assert configuration control.
	REQUIRE(true);
}
