// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file xml_parser_utils.h
 * @brief Core utility layer for parsing FERS XML scenario files.
 *
 * This file provides the internal mechanisms and data structures required to parse
 * individual XML elements into their corresponding simulation objects. It defines
 * a context-driven parsing approach, separating the extraction of XML data from
 * the global simulation state and managing external file dependencies through
 * function hooks.
 */

#pragma once

#include <filesystem>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/parameters.h"
#include "core/sim_id.h"
#include "radar/schedule_period.h"
#include "serial/libxml_wrapper.h"

// Forward declarations to minimize include dependencies
namespace antenna
{
	class Antenna;
}
namespace fers_signal
{
	class RadarSignal;
}
namespace timing
{
	class Timing;
}
namespace radar
{
	class Receiver;
	class Target;
	class Transmitter;
	class Platform;
}
namespace core
{
	class World;
}

namespace serial::xml_parser_utils
{
	/**
	 * @struct ReferenceLookup
	 * @brief Holds maps to resolve string names to internal SimId references during XML parsing.
	 *
	 * XML documents often cross-reference entities by name (e.g., a transmitter references
	 * an antenna by its string name). This struct provides the lookup tables needed to link
	 * these entities using their generated `SimId`.
	 */
	struct ReferenceLookup
	{
		const std::unordered_map<std::string, SimId>* waveforms; ///< Map of waveform names to IDs.
		const std::unordered_map<std::string, SimId>* antennas; ///< Map of antenna names to IDs.
		const std::unordered_map<std::string, SimId>* timings; ///< Map of timing object names to IDs.
	};

	/**
	 * @struct AssetLoaders
	 * @brief Container for functions that load external file-backed assets.
	 *
	 * Asset loading operations (such as reading waveforms, antenna patterns, and target
	 * RCS data from disk) are delegated to these `std::function` hooks. This allows
	 * the parser to flexibly resolve external file references during the XML parsing phase.
	 */
	struct AssetLoaders
	{
		/// Hook to load a pulsed waveform from an external file.
		std::function<std::unique_ptr<fers_signal::RadarSignal>(const std::string& name,
																const std::filesystem::path& pulse_path, RealType power,
																RealType carrierFreq, SimId id)>
			loadWaveform;

		/// Hook to load an antenna pattern defined in a legacy XML format.
		std::function<std::unique_ptr<antenna::Antenna>(const std::string& name, const std::string& filename, SimId id)>
			loadXmlAntenna;

		/// Hook to load an antenna pattern from an HDF5 file.
		std::function<std::unique_ptr<antenna::Antenna>(const std::string& name, const std::string& filename, SimId id)>
			loadH5Antenna;

		/// Hook to load a target's Radar Cross Section (RCS) from a file.
		std::function<std::unique_ptr<radar::Target>(radar::Platform* platform, const std::string& name,
													 const std::string& filename, unsigned seed, SimId id)>
			loadFileTarget;
	};

	/**
	 * @struct ParserContext
	 * @brief Encapsulates the state required during the XML parsing process.
	 *
	 * This context object holds the intermediate simulation parameters, a pointer to
	 * the simulation world, the base directory for resolving relative paths, the
	 * master random number generator, and the asset loading hooks. It is passed
	 * through the parsing functions to aggregate the scenario definition.
	 */
	struct ParserContext
	{
		params::Parameters parameters; ///< An isolated copy of the simulation parameters being built.
		core::World* world = nullptr; ///< Pointer to the World where parsed objects are inserted.
		std::filesystem::path base_dir; ///< The directory of the main XML file (used to resolve relative asset paths).
		std::mt19937* master_seeder = nullptr; ///< RNG used to generate independent seeds for simulated objects.
		AssetLoaders loaders; ///< The injected asset loaders for external files.
		std::unordered_map<SimId, std::shared_ptr<timing::Timing>>
			timing_instances; ///< Shared timing instances keyed by prototype ID.
	};

	/**
	 * @brief Extracts a floating-point (RealType) value from a named child element.
	 * @param element The parent XML element.
	 * @param elementName The name of the child element to extract text from.
	 * @return The parsed floating-point value.
	 * @throws XmlException if the child element is missing or empty.
	 */
	RealType get_child_real_type(const XmlElement& element, const std::string& elementName);

	/**
	 * @brief Extracts a boolean value from a named attribute.
	 * @param element The XML element containing the attribute.
	 * @param attributeName The name of the attribute.
	 * @param defaultVal The value to return if the attribute is missing or invalid.
	 * @return The parsed boolean value, or the default if the attribute is missing or invalid.
	 */
	bool get_attribute_bool(const XmlElement& element, const std::string& attributeName, bool defaultVal);

	/**
	 * @brief Generates a unique SimId based on the requested object type.
	 * @param owner The name/description of the object requesting the ID (used for logging).
	 * @param type The category/type of the object.
	 * @return A newly generated SimId.
	 */
	SimId assign_id_from_attribute(const std::string& owner, ObjectType type);

	/**
	 * @brief Resolves an XML string reference into an internal SimId.
	 * @param element The XML element containing the string reference attribute.
	 * @param attributeName The name of the attribute containing the reference string.
	 * @param owner A description of the object making the reference (used for error messages).
	 * @param name_map The lookup table mapping string names to SimIds.
	 * @return The resolved SimId.
	 * @throws XmlException if the reference cannot be resolved or is missing.
	 */
	SimId resolve_reference_id(const XmlElement& element, const std::string& attributeName, const std::string& owner,
							   const std::unordered_map<std::string, SimId>& name_map);

	/**
	 * @brief Parses a schedule (active periods) for a transmitter or receiver.
	 * @param parent The parent XML element that might contain a `<schedule>` block.
	 * @param parentName Name of the parent for error logging.
	 * @param isPulsed True if the owning object operates in pulsed mode (used for PRI validation).
	 * @param pri The pulse repetition interval, if applicable.
	 * @return A vector of parsed and validated `SchedulePeriod` objects.
	 */
	std::vector<radar::SchedulePeriod> parseSchedule(const XmlElement& parent, const std::string& parentName,
													 bool isPulsed, RealType pri = 0.0);

	/**
	 * @brief Parses the `<parameters>` block into the isolated context parameters.
	 * @param parameters The `<parameters>` XML element.
	 * @param params_out The `Parameters` struct to mutate with parsed values.
	 */
	void parseParameters(const XmlElement& parameters, params::Parameters& params_out);

	/**
	 * @brief Parses a `<waveform>` block and adds it to the World.
	 * @param waveform The `<waveform>` XML element.
	 * @param ctx The current parser context.
	 */
	void parseWaveform(const XmlElement& waveform, ParserContext& ctx);

	/**
	 * @brief Parses a `<timing>` block and adds the prototype timing to the World.
	 * @param timing The `<timing>` XML element.
	 * @param ctx The current parser context.
	 */
	void parseTiming(const XmlElement& timing, ParserContext& ctx);

	/**
	 * @brief Parses an `<antenna>` block and adds it to the World.
	 * @param antenna The `<antenna>` XML element.
	 * @param ctx The current parser context.
	 */
	void parseAntenna(const XmlElement& antenna, ParserContext& ctx);

	/**
	 * @brief Parses a `<motionpath>` block and attaches it to a Platform.
	 * @param motionPath The `<motionpath>` XML element.
	 * @param platform The platform to modify.
	 */
	void parseMotionPath(const XmlElement& motionPath, radar::Platform* platform);

	/**
	 * @brief Parses a `<rotationpath>` block and attaches it to a Platform.
	 * @param rotation The `<rotationpath>` XML element.
	 * @param platform The platform to modify.
	 */
	void parseRotationPath(const XmlElement& rotation, radar::Platform* platform, params::RotationAngleUnit unit);

	/**
	 * @brief Parses a `<fixedrotation>` block and attaches it to a Platform.
	 * @param rotation The `<fixedrotation>` XML element.
	 * @param platform The platform to modify.
	 */
	void parseFixedRotation(const XmlElement& rotation, radar::Platform* platform, params::RotationAngleUnit unit);

	/**
	 * @brief Parses a `<transmitter>` block, resolves its dependencies, and adds it to the World.
	 * @param transmitter The `<transmitter>` XML element.
	 * @param platform The platform this transmitter belongs to.
	 * @param ctx The current parser context.
	 * @param refs Lookup tables for resolving waveform, antenna, and timing references.
	 * @return A pointer to the newly created Transmitter object.
	 */
	radar::Transmitter* parseTransmitter(const XmlElement& transmitter, radar::Platform* platform, ParserContext& ctx,
										 const ReferenceLookup& refs);

	/**
	 * @brief Parses a `<receiver>` block, resolves its dependencies, and adds it to the World.
	 * @param receiver The `<receiver>` XML element.
	 * @param platform The platform this receiver belongs to.
	 * @param ctx The current parser context.
	 * @param refs Lookup tables for resolving antenna and timing references.
	 * @return A pointer to the newly created Receiver object.
	 */
	radar::Receiver* parseReceiver(const XmlElement& receiver, radar::Platform* platform, ParserContext& ctx,
								   const ReferenceLookup& refs);

	/**
	 * @brief Parses a `<monostatic>` block, creating a linked transmitter and receiver pair.
	 * @param monostatic The `<monostatic>` XML element.
	 * @param platform The platform this radar belongs to.
	 * @param ctx The current parser context.
	 * @param refs Lookup tables for resolving references.
	 */
	void parseMonostatic(const XmlElement& monostatic, radar::Platform* platform, ParserContext& ctx,
						 const ReferenceLookup& refs);

	/**
	 * @brief Parses a `<target>` block and adds it to the World.
	 * @param target The `<target>` XML element.
	 * @param platform The platform this target belongs to.
	 * @param ctx The current parser context.
	 */
	void parseTarget(const XmlElement& target, radar::Platform* platform, ParserContext& ctx);

	/**
	 * @brief Iterates and parses all children elements (radars, targets) of a platform.
	 * @param platform The `<platform>` XML element.
	 * @param ctx The current parser context.
	 * @param plat The Platform object to attach parsed elements to.
	 * @param register_name Callback used to ensure unique naming globally across parsed objects.
	 * @param refs Lookup tables for resolving references.
	 */
	void parsePlatformElements(const XmlElement& platform, ParserContext& ctx, radar::Platform* plat,
							   const std::function<void(const XmlElement&, std::string_view)>& register_name,
							   const ReferenceLookup& refs);

	/**
	 * @brief Parses a complete `<platform>` block, including its motion paths and sub-elements.
	 * @param platform The `<platform>` XML element.
	 * @param ctx The current parser context.
	 * @param register_name Callback used to ensure unique naming.
	 * @param refs Lookup tables for resolving references.
	 */
	void parsePlatform(const XmlElement& platform, ParserContext& ctx,
					   const std::function<void(const XmlElement&, std::string_view)>& register_name,
					   const ReferenceLookup& refs);

	/**
	 * @brief Recursively finds all `<include>` tags in a document and resolves their absolute paths.
	 * @param doc The XML document to search.
	 * @param currentDir The base directory for resolving relative paths.
	 * @param includePaths A vector populated with the absolute paths of included files.
	 */
	void collectIncludeElements(const XmlDocument& doc, const std::filesystem::path& currentDir,
								std::vector<std::filesystem::path>& includePaths);

	/**
	 * @brief Resolves and merges all `<include>` files directly into the provided main document.
	 * @param mainDoc The primary XML document that will be mutated.
	 * @param currentDir The base directory used to resolve include paths.
	 * @return True if at least one file was included and merged, false otherwise.
	 */
	bool addIncludeFilesToMainDocument(const XmlDocument& mainDoc, const std::filesystem::path& currentDir);

	/**
	 * @brief Validates an XML document against the embedded DTD and XSD schemas.
	 * @param didCombine Flag indicating whether the document contains merged includes (used for formatting log
	 * messages).
	 * @param mainDoc The XML document to validate.
	 * @throws XmlException if validation fails.
	 */
	void validateXml(bool didCombine, const XmlDocument& mainDoc);

	/**
	 * @brief Coordinates the full parsing of a validated XML document tree.
	 *
	 * This is the root parsing function that iterates over parameters, waveforms, timings,
	 * antennas, and platforms. It populates the World in the proper order and triggers
	 * initial event scheduling.
	 *
	 * @param doc The parsed XML document tree.
	 * @param ctx The parser context containing the World, isolated parameters, and asset loaders.
	 */
	void processParsedDocument(const XmlDocument& doc, ParserContext& ctx);

	/**
	 * @brief Creates an `AssetLoaders` struct populated with standard file-I/O implementations.
	 *
	 * Provides the default hooks to load actual HDF5, XML, and binary waveform files
	 * from the filesystem into the simulation environment.
	 *
	 * @return An `AssetLoaders` instance with standard file handlers attached.
	 */
	AssetLoaders createDefaultAssetLoaders();
}
