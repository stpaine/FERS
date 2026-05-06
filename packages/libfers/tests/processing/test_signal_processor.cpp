#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <memory>
#include <numeric>
#include <random>
#include <span>
#include <vector>

#include "core/config.h"
#include "core/parameters.h"
#include "interpolation/interpolation_point.h"
#include "processing/signal_processor.h"
#include "radar/platform.h"
#include "radar/transmitter.h"
#include "serial/response.h"
#include "signal/radar_signal.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};

	struct TestSignal final : public fers_signal::Signal
	{
		std::vector<ComplexType> data;

		std::vector<ComplexType> render(const std::vector<interp::InterpPoint>&, unsigned& size,
										RealType) const override
		{
			size = static_cast<unsigned>(data.size());
			return data;
		}
	};

	std::unique_ptr<serial::Response> makeResponse(const fers_signal::RadarSignal& wave,
												   const radar::Transmitter& transmitter,
												   const std::vector<interp::InterpPoint>& points)
	{
		auto response = std::make_unique<serial::Response>(&wave, &transmitter);
		for (const auto& point : points)
		{
			response->addInterpPoint(point);
		}
		return response;
	}

	RealType meanOfChannel(const std::vector<ComplexType>& window, bool realChannel)
	{
		RealType sum = 0.0;
		for (const auto& sample : window)
		{
			sum += realChannel ? sample.real() : sample.imag();
		}
		return sum / static_cast<RealType>(window.size());
	}

	RealType varianceOfChannel(const std::vector<ComplexType>& window, bool realChannel, RealType mean)
	{
		RealType accum = 0.0;
		for (const auto& sample : window)
		{
			const RealType value = realChannel ? sample.real() : sample.imag();
			const RealType diff = value - mean;
			accum += diff * diff;
		}
		return accum / static_cast<RealType>(window.size());
	}
}

TEST_CASE("renderWindow accumulates overlapping responses with offsets", "[processing][signal]")
{
	ParamGuard guard;
	params::setRate(4.0);
	params::setOversampleRatio(1);

	radar::Platform platform("RxPlatform");
	radar::Transmitter transmitter(&platform, "TxA", radar::OperationMode::PULSED_MODE, 11);

	auto signal_a = std::make_unique<TestSignal>();
	signal_a->data = {ComplexType{1.0, 0.0}, ComplexType{2.0, 1.0}, ComplexType{3.0, -1.0}};
	fers_signal::RadarSignal wave_a("wave_a", 1.0, 1.0, 1.0, std::move(signal_a), 101);

	auto signal_b = std::make_unique<TestSignal>();
	signal_b->data = {ComplexType{1.0, 1.0}, ComplexType{1.0, -1.0}};
	fers_signal::RadarSignal wave_b("wave_b", 1.0, 1.0, 1.0, std::move(signal_b), 102);

	auto signal_c = std::make_unique<TestSignal>();
	signal_c->data = {ComplexType{9.0, 9.0}};
	fers_signal::RadarSignal wave_c("wave_c", 1.0, 1.0, 1.0, std::move(signal_c), 103);

	std::vector<std::unique_ptr<serial::Response>> responses;
	responses.emplace_back(makeResponse(wave_a, transmitter, {{1.0, -0.25, 0.0, 0.0}, {1.0, 0.1, 0.0, 0.0}}));
	responses.emplace_back(makeResponse(wave_b, transmitter, {{1.0, 0.25, 0.0, 0.0}, {1.0, 0.5, 0.0, 0.0}}));
	responses.emplace_back(makeResponse(wave_c, transmitter, {{1.0, 2.0, 0.0, 0.0}, {1.0, 3.0, 0.0, 0.0}}));

	const RealType length = 1.0;
	const RealType start = 0.0;
	const RealType frac_delay = 0.0;
	const auto local_window_size = static_cast<unsigned>(std::ceil(length * params::rate()));
	std::vector<ComplexType> window(local_window_size, ComplexType{0.5, 0.0});

	processing::renderWindow(window, length, start, frac_delay, responses);

	REQUIRE(window.size() == 4);
	REQUIRE_THAT(window[0].real(), WithinAbs(0.5, 1e-12));
	REQUIRE_THAT(window[0].imag(), WithinAbs(0.0, 1e-12));

	REQUIRE_THAT(window[1].real(), WithinAbs(3.5, 1e-12));
	REQUIRE_THAT(window[1].imag(), WithinAbs(2.0, 1e-12));

	REQUIRE_THAT(window[2].real(), WithinAbs(4.5, 1e-12));
	REQUIRE_THAT(window[2].imag(), WithinAbs(-2.0, 1e-12));

	REQUIRE_THAT(window[3].real(), WithinAbs(0.5, 1e-12));
	REQUIRE_THAT(window[3].imag(), WithinAbs(0.0, 1e-12));
}

TEST_CASE("applyThermalNoise respects zero temperature", "[processing][signal]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);

	std::vector<ComplexType> window = {ComplexType{1.0, 2.0}, ComplexType{-3.0, 4.0}};
	std::mt19937 rng(1234u);

	processing::applyThermalNoise(window, 0.0, rng);

	REQUIRE_THAT(window[0].real(), WithinAbs(1.0, 0.0));
	REQUIRE_THAT(window[0].imag(), WithinAbs(2.0, 0.0));
	REQUIRE_THAT(window[1].real(), WithinAbs(-3.0, 0.0));
	REQUIRE_THAT(window[1].imag(), WithinAbs(4.0, 0.0));
}

TEST_CASE("applyThermalNoise produces expected power", "[processing][signal]")
{
	ParamGuard guard;
	params::setRate(2000.0);
	params::setOversampleRatio(2);

	const RealType temperature = 300.0;
	const RealType bandwidth = params::rate() / (2.0 * params::oversampleRatio());
	const RealType total_power = params::boltzmannK() * temperature * bandwidth;
	const RealType per_channel_power = total_power / 2.0;
	const RealType stddev = std::sqrt(per_channel_power);

	std::vector<ComplexType> window(50000, ComplexType{0.0, 0.0});
	std::mt19937 rng(42u);

	processing::applyThermalNoise(window, temperature, rng);

	const RealType mean_real = meanOfChannel(window, true);
	const RealType mean_imag = meanOfChannel(window, false);
	const RealType var_real = varianceOfChannel(window, true, mean_real);
	const RealType var_imag = varianceOfChannel(window, false, mean_imag);

	const RealType mean_tolerance = 5.0 * stddev / std::sqrt(static_cast<RealType>(window.size()));
	REQUIRE_THAT(mean_real, WithinAbs(0.0, mean_tolerance));
	REQUIRE_THAT(mean_imag, WithinAbs(0.0, mean_tolerance));
	REQUIRE_THAT(var_real, WithinRel(per_channel_power, 0.1));
	REQUIRE_THAT(var_imag, WithinRel(per_channel_power, 0.1));
}

TEST_CASE("applyThermalNoiseAtSampleRate scales complex variance to IF rate", "[processing][signal][fmcw][if]")
{
	ParamGuard guard;

	const RealType temperature = 290.0;
	const RealType if_rate = 64'000.0;
	const RealType complex_power = params::boltzmannK() * temperature * if_rate;
	const RealType per_channel_power = complex_power / 2.0;
	const RealType stddev = std::sqrt(per_channel_power);

	std::vector<ComplexType> window(60000, ComplexType{0.0, 0.0});
	std::mt19937 rng(123u);

	processing::applyThermalNoiseAtSampleRate(window, temperature, rng, if_rate);

	const RealType mean_real = meanOfChannel(window, true);
	const RealType mean_imag = meanOfChannel(window, false);
	const RealType var_real = varianceOfChannel(window, true, mean_real);
	const RealType var_imag = varianceOfChannel(window, false, mean_imag);

	const RealType mean_tolerance = 5.0 * stddev / std::sqrt(static_cast<RealType>(window.size()));
	REQUIRE_THAT(mean_real, WithinAbs(0.0, mean_tolerance));
	REQUIRE_THAT(mean_imag, WithinAbs(0.0, mean_tolerance));
	REQUIRE_THAT(var_real + var_imag, WithinRel(complex_power, 0.1));
}

TEST_CASE("quantizeAndScaleWindow normalizes when adc disabled", "[processing][signal]")
{
	ParamGuard guard;
	params::setAdcBits(0);

	std::vector<ComplexType> window = {ComplexType{2.0, 0.0}, ComplexType{-1.0, -3.0}};

	const RealType fullscale = processing::quantizeAndScaleWindow(window);

	REQUIRE_THAT(fullscale, WithinAbs(3.0, 1e-12));
	REQUIRE_THAT(window[0].real(), WithinAbs(2.0 / 3.0, 1e-12));
	REQUIRE_THAT(window[0].imag(), WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(window[1].real(), WithinAbs(-1.0 / 3.0, 1e-12));
	REQUIRE_THAT(window[1].imag(), WithinAbs(-1.0, 1e-12));
}

TEST_CASE("quantizeAndScaleWindow applies adc quantization", "[processing][signal]")
{
	ParamGuard guard;
	params::setAdcBits(3);

	std::vector<ComplexType> window = {ComplexType{2.0, -2.0}, ComplexType{1.0, -0.5}};

	const RealType fullscale = processing::quantizeAndScaleWindow(window);

	REQUIRE_THAT(fullscale, WithinAbs(2.0, 1e-12));
	REQUIRE_THAT(window[0].real(), WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(window[0].imag(), WithinAbs(-1.0, 1e-12));
	REQUIRE_THAT(window[1].real(), WithinAbs(0.5, 1e-12));
	REQUIRE_THAT(window[1].imag(), WithinAbs(-0.25, 1e-12));
}

TEST_CASE("quantizeAndScaleWindow leaves zeros unchanged", "[processing][signal]")
{
	ParamGuard guard;
	params::setAdcBits(0);

	std::vector<ComplexType> window = {ComplexType{0.0, 0.0}, ComplexType{0.0, 0.0}};

	const RealType fullscale = processing::quantizeAndScaleWindow(window);

	REQUIRE_THAT(fullscale, WithinAbs(0.0, 0.0));
	REQUIRE_THAT(window[0].real(), WithinAbs(0.0, 0.0));
	REQUIRE_THAT(window[0].imag(), WithinAbs(0.0, 0.0));
	REQUIRE_THAT(window[1].real(), WithinAbs(0.0, 0.0));
	REQUIRE_THAT(window[1].imag(), WithinAbs(0.0, 0.0));
}
