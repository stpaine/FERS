#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <libfers/api.h>
#include <nlohmann/json.hpp>
#include <string>

#include "api_test_helpers.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace
{
	struct PreviewIds
	{
		std::uint64_t tx = 0;
		std::uint64_t rx = 0;
		std::uint64_t target = 0;
	};

	template <std::size_t N>
	[[nodiscard]] std::string linkLabelString(const char (&label)[N])
	{
		const auto label_end = std::find(std::begin(label), std::end(label), '\0');
		return std::string(std::begin(label), label_end);
	}

	void requirePreviewIds(const PreviewIds& ids)
	{
		REQUIRE(ids.tx != 0);
		REQUIRE(ids.rx != 0);
		REQUIRE(ids.target != 0);
	}

	[[nodiscard]] PreviewIds previewIdsFromBistaticScenario(const nlohmann::json& scenario)
	{
		PreviewIds ids;
		for (const auto& platform : scenario.at("simulation").at("platforms"))
		{
			for (const auto& component : platform.at("components"))
			{
				if (component.contains("transmitter"))
				{
					ids.tx = api_test::parseId(component.at("transmitter").at("id"));
				}
				else if (component.contains("receiver"))
				{
					ids.rx = api_test::parseId(component.at("receiver").at("id"));
				}
				else if (component.contains("target"))
				{
					ids.target = api_test::parseId(component.at("target").at("id"));
				}
			}
		}
		return ids;
	}

	[[nodiscard]] PreviewIds previewIdsFromMonostaticScenario(const nlohmann::json& scenario)
	{
		PreviewIds ids;
		for (const auto& platform : scenario.at("simulation").at("platforms"))
		{
			for (const auto& component : platform.at("components"))
			{
				if (component.contains("monostatic"))
				{
					ids.tx = api_test::parseId(component.at("monostatic").at("tx_id"));
					ids.rx = api_test::parseId(component.at("monostatic").at("rx_id"));
				}
				else if (component.contains("target"))
				{
					ids.target = api_test::parseId(component.at("target").at("id"));
				}
			}
		}
		return ids;
	}

	void requireCommonPreviewLinkFields(const fers_visual_link_t& link)
	{
		const std::string label = linkLabelString(link.label);
		REQUIRE(link.label[sizeof(link.label) - 1] == '\0');
		REQUIRE_FALSE(label.empty());
		REQUIRE(std::isfinite(link.display_value));
	}

	void requireBistaticTxTargetLink(const fers_visual_link_t& link, const PreviewIds& ids)
	{
		const std::string label = linkLabelString(link.label);
		REQUIRE(link.source_id == ids.tx);
		REQUIRE(link.dest_id == ids.target);
		REQUIRE(link.origin_id == ids.tx);
		REQUIRE(link.quality == FERS_LINK_STRONG);
		REQUIRE_THAT(link.display_value, WithinAbs(std::stod(label), 0.1));
	}

	void requireDirectTxRxLink(const fers_visual_link_t& link, const PreviewIds& ids)
	{
		REQUIRE(link.source_id == ids.tx);
		REQUIRE(link.dest_id == ids.rx);
		REQUIRE(link.origin_id == ids.tx);
		REQUIRE(link.quality == FERS_LINK_STRONG);
	}

	void requireBistaticTargetRxLink(const fers_visual_link_t& link, const PreviewIds& ids)
	{
		const std::string label = linkLabelString(link.label);
		REQUIRE(link.source_id == ids.target);
		REQUIRE(link.dest_id == ids.rx);
		REQUIRE(link.origin_id == ids.tx);
		REQUIRE(link.quality == FERS_LINK_STRONG);
		REQUIRE_THAT(link.display_value, WithinAbs(std::stod(label), 0.1));
	}

	void requireBistaticPreviewLinks(const fers_visual_link_list_t& links, const PreviewIds& ids)
	{
		bool saw_tx_tgt = false;
		bool saw_direct = false;
		bool saw_tgt_rx = false;

		for (size_t i = 0; i < links.count; ++i)
		{
			const auto& link = links.links[i];
			requireCommonPreviewLinkFields(link);

			switch (link.type)
			{
			case FERS_LINK_BISTATIC_TX_TGT:
				saw_tx_tgt = true;
				requireBistaticTxTargetLink(link, ids);
				break;
			case FERS_LINK_DIRECT_TX_RX:
				saw_direct = true;
				requireDirectTxRxLink(link, ids);
				break;
			case FERS_LINK_BISTATIC_TGT_RX:
				saw_tgt_rx = true;
				requireBistaticTargetRxLink(link, ids);
				break;
			default:
				break;
			}
		}

		REQUIRE(saw_tx_tgt);
		REQUIRE(saw_direct);
		REQUIRE(saw_tgt_rx);
	}

	void requireMonostaticLink(const fers_visual_link_t& link, const PreviewIds& ids)
	{
		REQUIRE(link.source_id == ids.tx);
		REQUIRE(link.dest_id == ids.target);
		REQUIRE(link.origin_id == ids.tx);
		REQUIRE(link.quality == FERS_LINK_STRONG);
	}

	void requireMonostaticPreviewLinks(const fers_visual_link_list_t& links, const PreviewIds& ids)
	{
		bool saw_monostatic = false;
		for (size_t i = 0; i < links.count; ++i)
		{
			const auto& link = links.links[i];
			if (link.type == FERS_LINK_MONOSTATIC)
			{
				saw_monostatic = true;
				requireMonostaticLink(link, ids);
			}
		}

		REQUIRE(saw_monostatic);
	}
}

TEST_CASE("API antenna pattern validates context and antenna existence", "[api][preview]")
{
	api_test::clearLastError();

	SECTION("null context")
	{
		api_test::AntennaPattern const pattern(fers_get_antenna_pattern(nullptr, 1, 4, 3, 1.0e9));
		REQUIRE(pattern.get() == nullptr);
		api_test::ApiString const error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("context must be non-null"));
	}

	SECTION("unknown antenna id")
	{
		api_test::Context const context;
		REQUIRE(context.get() != nullptr);
		const std::string xml = api_test::previewScenarioXml();
		REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

		api_test::AntennaPattern const pattern(fers_get_antenna_pattern(context.get(), 999999, 4, 3, 1.0e9));
		REQUIRE(pattern.get() == nullptr);
		api_test::ApiString const error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("not found in the world"));
	}
}

TEST_CASE("API antenna pattern normalizes isotropic gains", "[api][preview]")
{
	api_test::ParamGuard const guard;
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::previewScenarioXml();
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	const auto scenario = api_test::parseScenarioJson(context.get());
	const auto antenna_id = api_test::parseId(scenario.at("simulation").at("antennas").at(0).at("id"));

	api_test::AntennaPattern const pattern(fers_get_antenna_pattern(context.get(), antenna_id, 4, 3, 1.0e9));
	REQUIRE(pattern.get() != nullptr);
	REQUIRE(pattern.get()->az_count == 4u);
	REQUIRE(pattern.get()->el_count == 3u);
	REQUIRE(pattern.get()->max_gain > 0.0);

	for (size_t i = 0; i < pattern.get()->az_count * pattern.get()->el_count; ++i)
	{
		REQUIRE_THAT(pattern.get()->gains[i], WithinAbs(1.0, 1e-12));
	}
}

TEST_CASE("API antenna pattern accepts zero frequency for isotropic previews", "[api][preview]")
{
	api_test::ParamGuard const guard;
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::previewScenarioXml();
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	const auto scenario = api_test::parseScenarioJson(context.get());
	const auto antenna_id = api_test::parseId(scenario.at("simulation").at("antennas").at(0).at("id"));

	api_test::AntennaPattern const pattern(fers_get_antenna_pattern(context.get(), antenna_id, 2, 2, 0.0));
	REQUIRE(pattern.get() != nullptr);
	REQUIRE(pattern.get()->max_gain > 0.0);
}

TEST_CASE("API antenna pattern previews XML antennas loaded from standalone files", "[api][preview]")
{
	api_test::ParamGuard const guard;
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	const api_test::ScopedPath antenna_path(api_test::uniqueTempPath("api_preview_xml_antenna", ".xml"));
	api_test::writeTextFile(antenna_path.path,
							R"(<?xml version="1.0" encoding="UTF-8"?>
<antenna>
  <azimuth unit="deg" format="dBi" symmetry="none">
    <gainsample><angle>-90</angle><gain>10</gain></gainsample>
    <gainsample><angle>0</angle><gain>0</gain></gainsample>
    <gainsample><angle>90</angle><gain>-10</gain></gainsample>
  </azimuth>
  <elevation unit="deg" format="dBi" symmetry="mirrored">
    <gainsample><angle>0</angle><gain>0</gain></gainsample>
    <gainsample><angle>90</angle><gain>0</gain></gainsample>
  </elevation>
</antenna>)");

	auto xml = api_test::previewScenarioXml("API XML Antenna Preview");
	const auto replace_all = [](std::string& text, const std::string_view from, const std::string& to)
	{
		size_t position = 0;
		while ((position = text.find(from, position)) != std::string::npos)
		{
			text.replace(position, from.size(), to);
			position += to.size();
		}
	};

	replace_all(xml, "<antenna name=\"api_preview_iso\" pattern=\"isotropic\"/>",
				"<antenna name=\"api_preview_xml\" pattern=\"xml\" filename=\"" + antenna_path.path.string() + "\"/>");
	replace_all(xml, "antenna=\"api_preview_iso\"", "antenna=\"api_preview_xml\"");

	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	const auto scenario = api_test::parseScenarioJson(context.get());
	const auto antenna_id = api_test::parseId(scenario.at("simulation").at("antennas").at(0).at("id"));

	api_test::AntennaPattern const pattern(fers_get_antenna_pattern(context.get(), antenna_id, 5, 3, 1.0e9));
	REQUIRE(pattern.get() != nullptr);
	REQUIRE(pattern.get()->max_gain > 0.0);
	REQUIRE_THAT(pattern.get()->gains[1 + 5], WithinAbs(1.0, 1e-12));
	REQUIRE(pattern.get()->gains[1 + 5] > pattern.get()->gains[3 + 5]);
}

TEST_CASE("API preview links return empty list for empty world", "[api][preview]")
{
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	api_test::PreviewLinks const links(fers_calculate_preview_links(context.get(), 0.0));
	REQUIRE(links.get() != nullptr);
	REQUIRE(links.get()->count == 0u);
	REQUIRE(links.get()->links == nullptr);
}

TEST_CASE("API preview links map link metadata into C structs", "[api][preview]")
{
	api_test::ParamGuard const guard;
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::previewScenarioXml();
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	const auto scenario = api_test::parseScenarioJson(context.get());
	const PreviewIds ids = previewIdsFromBistaticScenario(scenario);
	requirePreviewIds(ids);

	api_test::PreviewLinks const links(fers_calculate_preview_links(context.get(), 0.0));
	REQUIRE(links.get() != nullptr);
	REQUIRE(links.get()->count > 0u);
	requireBistaticPreviewLinks(*links.get(), ids);
}

TEST_CASE("API preview links map monostatic link types into C enums", "[api][preview]")
{
	api_test::ParamGuard const guard;
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::monostaticPreviewScenarioXml();
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	const auto scenario = api_test::parseScenarioJson(context.get());
	const PreviewIds ids = previewIdsFromMonostaticScenario(scenario);
	requirePreviewIds(ids);

	api_test::PreviewLinks const links(fers_calculate_preview_links(context.get(), 0.0));
	REQUIRE(links.get() != nullptr);
	requireMonostaticPreviewLinks(*links.get(), ids);
}
