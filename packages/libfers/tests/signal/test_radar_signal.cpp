#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <optional>
#include <span>
#include <vector>

#include "core/parameters.h"
#include "core/sim_id.h"
#include "interpolation/interpolation_filter.h"
#include "interpolation/interpolation_point.h"
#include "signal/radar_signal.h"

using Catch::Matchers::WithinAbs;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};

	ComplexType expectedConstantRender(const std::span<const RealType> filter, const RealType amplitude,
									   const RealType phase)
	{
		const int filt_length = static_cast<int>(filter.size());
		RealType sum = 0.0;
		for (int j = -filt_length / 2; j < filt_length / 2; ++j)
		{
			const int idx = j + filt_length / 2;
			if (idx >= 0 && idx < filt_length)
			{
				sum += filter[static_cast<size_t>(idx)];
			}
		}
		return amplitude * std::exp(ComplexType{0.0, 1.0} * phase) * sum;
	}
}

TEST_CASE("CwSignal render returns empty data", "[signal][radar][cw]")
{
	fers_signal::CwSignal signal;
	unsigned size = 99;
	const std::vector<interp::InterpPoint> points = {{1.0, 0.0, 0.0, 0.0}};

	auto data = signal.render(points, size, 0.0);

	REQUIRE(data.empty());
	REQUIRE(size == 0);
}

TEST_CASE("RadarSignal requires a signal", "[signal][radar]")
{
	REQUIRE_THROWS_AS(fers_signal::RadarSignal("test", 1.0, 1.0, 1.0, std::unique_ptr<fers_signal::Signal>{}, 0),
					  std::runtime_error);
}

TEST_CASE("FmcwChirpSignal applies sweep direction to phase only", "[signal][radar][fmcw]")
{
	const RealType bandwidth = 2.0e6;
	const RealType duration = 1.0e-3;
	const RealType offset = 1.0e5;
	const RealType u = 4.0e-4;

	fers_signal::FmcwChirpSignal up(bandwidth, duration, duration, offset);
	fers_signal::FmcwChirpSignal down(bandwidth, duration, duration, offset, std::nullopt,
									  fers_signal::FmcwChirpDirection::Down);

	const RealType alpha = bandwidth / duration;
	REQUIRE_FALSE(up.isDownChirp());
	REQUIRE(down.isDownChirp());
	REQUIRE_THAT(up.getChirpRate(), WithinAbs(alpha, 1e-6));
	REQUIRE_THAT(down.getChirpRate(), WithinAbs(alpha, 1e-6));
	REQUIRE_THAT(up.getSignedChirpRate(), WithinAbs(alpha, 1e-6));
	REQUIRE_THAT(down.getSignedChirpRate(), WithinAbs(-alpha, 1e-6));
	REQUIRE_THAT(up.basebandPhaseForChirpTime(u), WithinAbs(2.0 * PI * offset * u + PI * alpha * u * u, 1e-9));
	REQUIRE_THAT(down.basebandPhaseForChirpTime(u), WithinAbs(2.0 * PI * offset * u - PI * alpha * u * u, 1e-9));
}

TEST_CASE("FmcwTriangleSignal keeps phase continuous at leg and period boundaries", "[signal][radar][fmcw]")
{
	const RealType bandwidth = 2.0e6;
	const RealType duration = 1.0e-3;
	const RealType offset = 1.0e5;
	const RealType alpha = bandwidth / duration;
	const RealType eps = 1.0e-10;

	fers_signal::FmcwTriangleSignal triangle(bandwidth, duration, offset, 4);

	REQUIRE(triangle.isFmcwFamily());
	REQUIRE(triangle.isTriangle());
	REQUIRE_THAT(triangle.getChirpRate(), WithinAbs(alpha, 1e-6));
	REQUIRE_THAT(triangle.getTrianglePeriod(), WithinAbs(2.0 * duration, 1e-15));
	REQUIRE_THAT(triangle.getDeltaPhiUp(),
				 WithinAbs(2.0 * PI * offset * duration + PI * alpha * duration * duration, 1e-6));

	const RealType apex_left = triangle.basebandPhaseForTriangleTime(duration - eps);
	const RealType apex = triangle.basebandPhaseForTriangleTime(duration);
	const RealType apex_right = triangle.basebandPhaseForTriangleTime(duration + eps);
	REQUIRE_THAT(apex, WithinAbs(triangle.getDeltaPhiUp(), 1e-6));
	REQUIRE(std::abs(apex - apex_left) < 2.0);
	REQUIRE(std::abs(apex_right - apex) < 2.0);

	const RealType period_phase = triangle.basebandPhaseForTriangleTime(triangle.getTrianglePeriod());
	REQUIRE_THAT(period_phase, WithinAbs(2.0 * triangle.getDeltaPhiUp(), 1e-6));
	REQUIRE(triangle.instantaneousBasebandPhase(3.5 * triangle.getTrianglePeriod()).has_value());
	REQUIRE_FALSE(triangle.instantaneousBasebandPhase(4.0 * triangle.getTrianglePeriod()).has_value());
}

TEST_CASE("RadarSignal exposes metadata", "[signal][radar]")
{
	ParamGuard guard;
	params::setOversampleRatio(1);
	std::vector<ComplexType> data = {ComplexType{1.0, 0.0}};
	auto signal = std::make_unique<fers_signal::Signal>();
	signal->load(data, static_cast<unsigned>(data.size()), 1000.0);

	fers_signal::RadarSignal radar("waveform", 9.0, 77.0, 0.01, std::move(signal), 42);
	radar.setFilename("waveform.bin");

	REQUIRE(radar.getName() == "waveform");
	REQUIRE(radar.getCarrier() == 77.0);
	REQUIRE(radar.getPower() == 9.0);
	REQUIRE(radar.getLength() == 0.01);
	REQUIRE(radar.getRate() == 1000.0);
	REQUIRE(radar.getId() == 42);
	REQUIRE(radar.getFilename().has_value());
	REQUIRE(radar.getFilename().value() == "waveform.bin");
}

TEST_CASE("RadarSignal autogenerates waveform ids", "[signal][radar]")
{
	auto signal = std::make_unique<fers_signal::Signal>();
	fers_signal::RadarSignal radar("waveform", 1.0, 1.0, 1.0, std::move(signal), 0);

	REQUIRE(SimIdGenerator::getType(radar.getId()) == ObjectType::Waveform);
}

TEST_CASE("RadarSignal scales rendered data by power", "[signal][radar]")
{
	struct TestSignal final : public fers_signal::Signal
	{
		std::vector<ComplexType> stored;
		std::vector<ComplexType> render(const std::vector<interp::InterpPoint>&, unsigned& size,
										RealType) const override
		{
			size = static_cast<unsigned>(stored.size());
			return stored;
		}
	};

	auto signal = std::make_unique<TestSignal>();
	signal->stored = {ComplexType{1.0, 0.0}, ComplexType{2.0, -1.0}};

	fers_signal::RadarSignal radar("scaled", 4.0, 1.0, 1.0, std::move(signal), 7);

	unsigned size = 0;
	const std::vector<interp::InterpPoint> points = {{1.0, 0.0, 0.0, 0.0}};
	auto data = radar.render(points, size, 0.0);

	REQUIRE(size == 2);
	REQUIRE_THAT(data[0].real(), WithinAbs(2.0, 1e-12));
	REQUIRE_THAT(data[0].imag(), WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(data[1].real(), WithinAbs(4.0, 1e-12));
	REQUIRE_THAT(data[1].imag(), WithinAbs(-2.0, 1e-12));
}

TEST_CASE("Signal load updates size and rate", "[signal][radar]")
{
	ParamGuard guard;
	params::setOversampleRatio(2);
	std::vector<ComplexType> input = {ComplexType{1.0, 0.0}, ComplexType{0.5, -0.5}};

	fers_signal::Signal signal;
	signal.load(input, static_cast<unsigned>(input.size()), 100.0);

	REQUIRE(signal.getRate() == 200.0);

	const std::vector<interp::InterpPoint> points = {{1.0, 0.0, 0.0, 0.0}};
	unsigned size = 0;
	const auto data = signal.render(points, size, 0.0);
	REQUIRE(size == input.size() * 2);
	REQUIRE(data.size() == size);
}

TEST_CASE("Signal render matches constant input physics", "[signal][radar]")
{
	ParamGuard guard;
	params::setOversampleRatio(1);

	const unsigned filter_length = params::renderFilterLength();
	const unsigned sample_count = filter_length * 4;
	std::vector<ComplexType> input(sample_count, ComplexType{1.0, 0.0});

	fers_signal::Signal signal;
	signal.load(input, sample_count, 1.0);

	const RealType power = 4.0;
	const RealType phase = PI / 4.0;
	const std::vector<interp::InterpPoint> points = {{power, 0.0, 0.0, phase}};

	unsigned size = 0;
	const auto data = signal.render(points, size, 0.0);

	const auto& interp = interp::InterpFilter::getInstance();
	const auto filter = interp.getFilter(0.0);
	const auto expected = expectedConstantRender(filter, std::sqrt(power), phase);

	const unsigned sample_index = filter_length;
	REQUIRE(sample_index < data.size());
	REQUIRE_THAT(data[sample_index].real(), WithinAbs(expected.real(), 1e-6));
	REQUIRE_THAT(data[sample_index].imag(), WithinAbs(expected.imag(), 1e-6));
}

TEST_CASE("Signal render interpolates power and phase", "[signal][radar]")
{
	ParamGuard guard;
	params::setOversampleRatio(1);

	const unsigned filter_length = params::renderFilterLength();
	const unsigned sample_count = filter_length * 4;
	std::vector<ComplexType> input(sample_count, ComplexType{1.0, 0.0});

	fers_signal::Signal signal;
	signal.load(input, sample_count, 1.0);

	const RealType power_a = 1.0;
	const RealType power_b = 9.0;
	const RealType phase_a = 0.0;
	const RealType phase_b = PI / 2.0;
	const std::vector<interp::InterpPoint> points = {{power_a, 0.0, 0.0, phase_a},
													 {power_b, 2.0 * filter_length, 0.0, phase_b}};

	unsigned size = 0;
	const auto data = signal.render(points, size, 0.0);

	const auto& interp = interp::InterpFilter::getInstance();
	const auto filter = interp.getFilter(0.0);

	const RealType bw = 0.5;
	const RealType amplitude = std::lerp(std::sqrt(power_a), std::sqrt(power_b), bw);
	const RealType phase = std::lerp(phase_a, phase_b, bw);
	const auto expected = expectedConstantRender(filter, amplitude, phase);

	const unsigned sample_index = filter_length;
	REQUIRE(sample_index < data.size());
	REQUIRE_THAT(data[sample_index].real(), WithinAbs(expected.real(), 1e-6));
	REQUIRE_THAT(data[sample_index].imag(), WithinAbs(expected.imag(), 1e-6));
}

TEST_CASE("Signal render responds to fractional delay", "[signal][radar]")
{
	ParamGuard guard;
	params::setOversampleRatio(1);

	const unsigned filter_length = params::renderFilterLength();
	const unsigned sample_count = filter_length * 4;
	std::vector<ComplexType> input(sample_count, ComplexType{1.0, 0.0});

	fers_signal::Signal signal;
	signal.load(input, sample_count, 1.0);

	const std::vector<interp::InterpPoint> points = {{1.0, 0.0, 0.0, 0.0}};

	unsigned size = 0;
	const auto data = signal.render(points, size, 0.25);

	const auto& interp = interp::InterpFilter::getInstance();
	const auto filter = interp.getFilter(0.75);
	const auto expected = expectedConstantRender(filter, 1.0, 0.0);

	const unsigned sample_index = filter_length;
	REQUIRE(sample_index < data.size());
	REQUIRE_THAT(data[sample_index].real(), WithinAbs(expected.real(), 1e-6));
	REQUIRE_THAT(data[sample_index].imag(), WithinAbs(expected.imag(), 1e-6));
}
