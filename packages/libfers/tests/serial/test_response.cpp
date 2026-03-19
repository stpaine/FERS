#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <memory>
#include <vector>

#include "core/config.h"
#include "interpolation/interpolation_point.h"
#include "radar/platform.h"
#include "radar/transmitter.h"
#include "serial/response.h"
#include "signal/radar_signal.h"

using Catch::Matchers::WithinAbs;

namespace
{
	struct CaptureSignal final : public fers_signal::Signal
	{
		std::vector<ComplexType> data;
		mutable std::vector<interp::InterpPoint> last_points;
		mutable RealType last_frac_delay = 0.0;

		std::vector<ComplexType> render(const std::vector<interp::InterpPoint>& points, unsigned& size,
										RealType fracWinDelay) const override
		{
			last_points = points;
			last_frac_delay = fracWinDelay;
			size = static_cast<unsigned>(data.size());
			return data;
		}
	};
}

TEST_CASE("Response start/end time default to zero", "[serial][response]")
{
	radar::Platform platform("TxPlatform");
	radar::Transmitter transmitter(&platform, "TxA", radar::OperationMode::PULSED_MODE, 7);

	auto signal = std::make_unique<CaptureSignal>();
	fers_signal::RadarSignal wave("wave", 1.0, 1.0, 1.0, std::move(signal), 100);

	serial::Response response(&wave, &transmitter);

	REQUIRE_THAT(response.startTime(), WithinAbs(0.0, 0.0));
	REQUIRE_THAT(response.endTime(), WithinAbs(0.0, 0.0));
	REQUIRE_THAT(response.getLength(), WithinAbs(0.0, 0.0));
}

TEST_CASE("Response records interpolation points", "[serial][response]")
{
	radar::Platform platform("TxPlatform");
	radar::Transmitter transmitter(&platform, "TxA", radar::OperationMode::PULSED_MODE, 7);

	auto signal = std::make_unique<CaptureSignal>();
	fers_signal::RadarSignal wave("wave", 1.0, 1.0, 1.0, std::move(signal), 100);

	serial::Response response(&wave, &transmitter);
	response.addInterpPoint({1.0, 0.25, 0.0, 0.0});
	response.addInterpPoint({4.0, 1.75, 0.0, 0.0});

	REQUIRE_THAT(response.startTime(), WithinAbs(0.25, 1e-12));
	REQUIRE_THAT(response.endTime(), WithinAbs(1.75, 1e-12));
	REQUIRE_THAT(response.getLength(), WithinAbs(1.5, 1e-12));
}

TEST_CASE("Response exposes transmitter id", "[serial][response]")
{
	radar::Platform platform("TxPlatform");
	radar::Transmitter transmitter(&platform, "TxA", radar::OperationMode::CW_MODE, 1234);

	auto signal = std::make_unique<CaptureSignal>();
	fers_signal::RadarSignal wave("wave", 1.0, 1.0, 1.0, std::move(signal), 100);

	serial::Response response(&wave, &transmitter);
	REQUIRE(response.getTransmitterId() == 1234);
}

TEST_CASE("Response renderBinary delegates to signal", "[serial][response]")
{
	radar::Platform platform("TxPlatform");
	radar::Transmitter transmitter(&platform, "TxA", radar::OperationMode::PULSED_MODE, 7);

	auto signal = std::make_unique<CaptureSignal>();
	signal->data = {ComplexType{1.0, -1.0}, ComplexType{0.5, 0.25}};
	std::vector<ComplexType> dummy = {ComplexType{0.0, 0.0}};
	signal->load(dummy, static_cast<unsigned>(dummy.size()), 1250.0);

	fers_signal::RadarSignal wave("wave", 1.0, 1.0, 1.0, std::move(signal), 100);

	serial::Response response(&wave, &transmitter);
	response.addInterpPoint({1.0, 0.0, 0.0, 0.0});
	response.addInterpPoint({2.0, 1.0, 0.5, 0.25});

	RealType rate = 0.0;
	unsigned size = 0;
	const RealType frac_delay = 0.125;
	const auto data = response.renderBinary(rate, size, frac_delay);

	REQUIRE_THAT(rate, WithinAbs(1250.0, 1e-12));
	REQUIRE(size == 2u);
	REQUIRE(data.size() == 2u);
	REQUIRE_THAT(data[0].real(), WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(data[0].imag(), WithinAbs(-1.0, 1e-12));
	REQUIRE_THAT(data[1].real(), WithinAbs(0.5, 1e-12));
	REQUIRE_THAT(data[1].imag(), WithinAbs(0.25, 1e-12));

	const auto* capture = dynamic_cast<const CaptureSignal*>(wave.getSignal());
	REQUIRE(capture != nullptr);
	REQUIRE(capture->last_points.size() == 2u);
	REQUIRE_THAT(capture->last_points[0].time, WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(capture->last_points[1].delay, WithinAbs(0.5, 1e-12));
	REQUIRE_THAT(capture->last_frac_delay, WithinAbs(frac_delay, 1e-12));
}
