/**
 * @file fers_context.h
 * @brief Internal C++ class that encapsulates the state of a simulation instance.
 */

#pragma once

#include <memory>
#include <random>

#include "world.h"

/**
 * @class FersContext
 * @brief Manages the lifetime and state of a single FERS simulation scenario.
 *
 * This class serves as the C++ backend for the opaque `fers_context_t`
 * handle exposed by the C-API. The primary reason for its existence is to
 * apply the "Pimpl" (Pointer to Implementation) idiom at the ABI boundary. By
 * hiding the C++ standard library types (`std::unique_ptr`, `std::mt19937`) and
 * the full definition of `core::World` from the C header, we create a stable
 * ABI that does not break when the internal C++ implementation changes. This
 * encapsulation ensures that clients of the C-API (like Rust) do not need to be
 * recompiled if only the library's internals are modified.
 *
 * Its secondary role is to own the `core::World` object, which represents the
 * entire scenario, and the master random number generator. This guarantees that
 * the simulation state persists in memory between API calls and that randomness
 * can be controlled deterministically from a single source.
 */
class FersContext
{
public:
	/**
	 * @brief Constructs a new simulation context, initializing an empty world.
	 *
	 * The master random number generator is default-constructed. This is a
	 * deliberate choice to allow the seed to be configured later, typically
	 * after parsing a scenario file. This ensures that the scenario itself can
	 * define its own seed for reproducible simulations. If no seed is provided,
	 * a random one is generated at load time.
	 */
	// NOLINTNEXTLINE(cert-msc51-cpp)
	FersContext() : _world(std::make_unique<core::World>()) {}

	/**
	 * @brief Retrieves a pointer to the simulation world.
	 * This provides direct, mutable access to the in-memory representation of the
	 * scenario, allowing API functions to modify it.
	 * @return A non-owning pointer to the `core::World` object.
	 */
	[[nodiscard]] core::World* getWorld() const noexcept { return _world.get(); }

	/**
	 * @brief Retrieves a mutable reference to the master random number seeder.
	 *
	 * A single master generator is used to seed all other random number
	 * generators within the simulation (e.g., for noise models, RCS
	 * fluctuations). This design is crucial for ensuring that a simulation can
	 * be made fully deterministic and reproducible by controlling a single seed
	 * value at the top level.
	 * @return A reference to the `std::mt19937` engine.
	 */
	[[nodiscard]] std::mt19937& getMasterSeeder() noexcept { return _master_seeder; }

private:
	/// Owns the `core::World` object, which contains all simulation entities.
	/// Using `std::unique_ptr` ensures that the world's complex state is
	/// automatically cleaned up when the FersContext is destroyed.
	std::unique_ptr<core::World> _world;

	/// Master random engine used to seed all other random generators in the simulation.
	std::mt19937 _master_seeder;
};
