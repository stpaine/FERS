// SPDX-License-Identifier: GPL-2.0-only
//
// Standalone oversampling audit harness for isolated script-style analysis.

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "core/logging.h"
#include "core/parameters.h"
#include "interpolation/interpolation_filter.h"
#include "interpolation/interpolation_point.h"
#include "processing/finalizer_pipeline.h"
#include "processing/signal_processor.h"
#include "serial/response.h"
#include "signal/dsp_filters.h"
#include "signal/radar_signal.h"

namespace fs = std::filesystem;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};

	struct SourceMapping
	{
		std::string suite;
		std::string upstream_file;
		std::string upstream_lines;
		std::string local_file;
		std::string note;
	};

	struct SuiteResult
	{
		std::string suite;
		fs::path output_dir;
		std::vector<std::string> generated_files;
		std::vector<std::pair<std::string, std::string>> summary_rows;
		std::vector<SourceMapping> source_mappings;
		std::string readme;
	};

	enum class CommandType
	{
		Suite,
		All,
		Matrix
	};

	struct CliOptions
	{
		CommandType command = CommandType::All;
		std::string suite = "all";
		unsigned ratio = 4;
		unsigned filter_length = 33;
		unsigned seed = 42;
		std::vector<unsigned> ratios = {1, 2, 4, 8};
		std::vector<unsigned> filter_lengths = {8, 33};
		fs::path output_dir = fs::path("verification/oversampling/output");
		fs::path executable_path;
	};

	struct RenderComparison
	{
		RealType rms_error{};
		RealType normalized_correlation{};
		RealType fitted_gain{};
		RealType fitted_phase{};
	};

	std::string toString(const RealType value, const int precision = 12)
	{
		std::ostringstream stream;
		stream << std::setprecision(precision) << std::scientific << value;
		return stream.str();
	}

	std::string toString(const unsigned value)
	{
		return std::to_string(value);
	}

	std::string jsonEscape(const std::string& value)
	{
		std::string out;
		out.reserve(value.size() + 8);
		for (const char ch : value)
		{
			switch (ch)
			{
			case '\\':
				out += "\\\\";
				break;
			case '"':
				out += "\\\"";
				break;
			case '\n':
				out += "\\n";
				break;
			case '\r':
				out += "\\r";
				break;
			case '\t':
				out += "\\t";
				break;
			default:
				out += ch;
				break;
			}
		}
		return out;
	}

	void ensureDirectory(const fs::path& directory)
	{
		std::error_code error;
		fs::create_directories(directory, error);
		if (error)
		{
			throw std::runtime_error("Failed to create directory: " + directory.string());
		}
	}

	void writeTextFile(const fs::path& path, const std::string& text)
	{
		ensureDirectory(path.parent_path());
		std::ofstream file(path);
		if (!file.is_open())
		{
			throw std::runtime_error("Failed to open file for writing: " + path.string());
		}
		file << text;
	}

	void writeCsv(const fs::path& path, const std::vector<std::string>& header,
				  const std::vector<std::vector<std::string>>& rows)
	{
		ensureDirectory(path.parent_path());
		std::ofstream file(path);
		if (!file.is_open())
		{
			throw std::runtime_error("Failed to open file for writing: " + path.string());
		}

		for (size_t i = 0; i < header.size(); ++i)
		{
			file << header[i];
			if (i + 1 < header.size())
			{
				file << ',';
			}
		}
		file << '\n';

		for (const auto& row : rows)
		{
			for (size_t i = 0; i < row.size(); ++i)
			{
				file << row[i];
				if (i + 1 < row.size())
				{
					file << ',';
				}
			}
			file << '\n';
		}
	}

	void addGeneratedFile(SuiteResult& result, const fs::path& path)
	{
		result.generated_files.push_back(fs::relative(path, result.output_dir).string());
	}

	std::string shellQuote(const std::string& value)
	{
		std::string out = "'";
		for (const char ch : value)
		{
			if (ch == '\'')
			{
				out += "'\"'\"'";
			}
			else
			{
				out += ch;
			}
		}
		out += "'";
		return out;
	}

	std::vector<unsigned> parseUnsignedList(const std::string& value)
	{
		std::vector<unsigned> out;
		std::stringstream stream(value);
		std::string token;
		while (std::getline(stream, token, ','))
		{
			if (!token.empty())
			{
				out.push_back(static_cast<unsigned>(std::stoul(token)));
			}
		}
		if (out.empty())
		{
			throw std::runtime_error("Expected a non-empty comma-separated unsigned list");
		}
		return out;
	}

	CliOptions parseCli(const int argc, char* argv[])
	{
		CliOptions options;
		options.executable_path = fs::absolute(argv[0]);

		if (argc >= 2)
		{
			const std::string command = argv[1];
			if (command == "suite")
			{
				options.command = CommandType::Suite;
				if (argc < 3)
				{
					throw std::runtime_error("suite command requires a suite name");
				}
				options.suite = argv[2];
			}
			else if (command == "matrix")
			{
				options.command = CommandType::Matrix;
			}
			else if (command == "all")
			{
				options.command = CommandType::All;
			}
			else
			{
				throw std::runtime_error("Unknown command: " + command);
			}
		}

		const int start_index = options.command == CommandType::Suite ? 3 : 2;
		for (int i = start_index; i < argc; ++i)
		{
			const std::string arg = argv[i];
			if (arg == "--ratio" && i + 1 < argc)
			{
				options.ratio = static_cast<unsigned>(std::stoul(argv[++i]));
			}
			else if (arg == "--filter-length" && i + 1 < argc)
			{
				options.filter_length = static_cast<unsigned>(std::stoul(argv[++i]));
			}
			else if (arg == "--seed" && i + 1 < argc)
			{
				options.seed = static_cast<unsigned>(std::stoul(argv[++i]));
			}
			else if (arg == "--output-dir" && i + 1 < argc)
			{
				options.output_dir = argv[++i];
			}
			else if (arg == "--ratios" && i + 1 < argc)
			{
				options.ratios = parseUnsignedList(argv[++i]);
			}
			else if (arg == "--filter-lengths" && i + 1 < argc)
			{
				options.filter_lengths = parseUnsignedList(argv[++i]);
			}
			else
			{
				throw std::runtime_error("Unknown or incomplete option: " + arg);
			}
		}

		return options;
	}

	std::vector<RealType> blackmanFirAudit(const RealType cutoff, const unsigned render_filter_length)
	{
		const unsigned filt_length = render_filter_length * 2;
		std::vector<RealType> coeffs(filt_length);
		const RealType n = filt_length / 2.0;
		const RealType pi_n = PI / n;

		for (unsigned i = 0; i < filt_length; ++i)
		{
			const RealType sinc_arg = cutoff * (static_cast<RealType>(i) - n);
			const RealType sinc_val = sinc_arg == 0 ? 1.0 : std::sin(sinc_arg * PI) / (sinc_arg * PI);
			const RealType window =
				0.42 - 0.5 * std::cos(pi_n * static_cast<RealType>(i)) + 0.08 * std::cos(2.0 * pi_n * static_cast<RealType>(i));
			coeffs[i] = sinc_val * window;
		}

		return coeffs;
	}

	std::vector<ComplexType> makeTone(const size_t sample_count, const RealType frequency)
	{
		std::vector<ComplexType> out(sample_count);
		for (size_t i = 0; i < sample_count; ++i)
		{
			const RealType phase = 2.0 * PI * frequency * static_cast<RealType>(i);
			out[i] = {std::cos(phase), std::sin(phase)};
		}
		return out;
	}

	std::vector<ComplexType> makeChirp(const size_t sample_count, const RealType start_frequency, const RealType end_frequency)
	{
		std::vector<ComplexType> out(sample_count);
		for (size_t i = 0; i < sample_count; ++i)
		{
			const RealType t = sample_count > 1 ? static_cast<RealType>(i) / static_cast<RealType>(sample_count - 1) : 0.0;
			const RealType frequency = std::lerp(start_frequency, end_frequency, t);
			const RealType phase = 2.0 * PI * frequency * static_cast<RealType>(i);
			out[i] = {std::cos(phase), std::sin(phase)};
		}
		return out;
	}

	std::vector<ComplexType> makeLowBandNoise(const size_t sample_count, std::mt19937& rng)
	{
		std::uniform_real_distribution<RealType> phase_dist(0.0, 2.0 * PI);
		std::uniform_real_distribution<RealType> freq_dist(0.01, 0.12);
		std::vector<RealType> freqs(16);
		std::vector<RealType> phases(16);
		for (size_t i = 0; i < freqs.size(); ++i)
		{
			freqs[i] = freq_dist(rng);
			phases[i] = phase_dist(rng);
		}

		std::vector<ComplexType> out(sample_count, ComplexType{0.0, 0.0});
		for (size_t n = 0; n < sample_count; ++n)
		{
			ComplexType accum{0.0, 0.0};
			for (size_t i = 0; i < freqs.size(); ++i)
			{
				const RealType phase = 2.0 * PI * freqs[i] * static_cast<RealType>(n) + phases[i];
				accum += ComplexType(std::cos(phase), std::sin(phase));
			}
			out[n] = accum / static_cast<RealType>(freqs.size());
		}
		return out;
	}

	RenderComparison compareSignals(std::span<const ComplexType> reference, std::span<const ComplexType> observed,
									const size_t skip_edge)
	{
		const size_t start = std::min(skip_edge, std::min(reference.size(), observed.size()));
		const size_t end = reference.size() > skip_edge && observed.size() > skip_edge ?
			std::min(reference.size(), observed.size()) - skip_edge :
			std::min(reference.size(), observed.size());

		ComplexType correlation{0.0, 0.0};
		RealType output_energy = 0.0;
		RealType reference_energy = 0.0;
		RealType mse = 0.0;
		size_t count = 0;

		for (size_t i = start; i < end; ++i)
		{
			correlation += observed[i] * std::conj(reference[i]);
			output_energy += std::norm(observed[i]);
			reference_energy += std::norm(reference[i]);
			mse += std::norm(observed[i] - reference[i]);
			++count;
		}

		RenderComparison result;
		result.rms_error = count == 0 ? 0.0 : std::sqrt(mse / static_cast<RealType>(count));
		result.normalized_correlation =
			(output_energy == 0.0 || reference_energy == 0.0) ? 0.0 : std::abs(correlation) / std::sqrt(output_energy * reference_energy);

		ComplexType fitted = reference_energy == 0.0 ? ComplexType{0.0, 0.0} : correlation / reference_energy;
		result.fitted_gain = std::abs(fitted);
		result.fitted_phase = std::arg(fitted);
		return result;
	}

	RealType filterTapSum(const std::span<const RealType> filter)
	{
		return std::accumulate(filter.begin(), filter.end(), 0.0);
	}

	ComplexType expectedConstantRender(const std::span<const RealType> filter, const RealType amplitude,
									   const RealType phase, const RealType upsample_dc_gain)
	{
		return amplitude * upsample_dc_gain * std::exp(ComplexType{0.0, 1.0} * phase) * filterTapSum(filter);
	}

	RealType wrappedDelayFromFracWinDelay(const RealType frac_win_delay)
	{
		RealType fdelay = -frac_win_delay;
		const int unwrap = static_cast<int>(std::floor(fdelay));
		fdelay -= unwrap;
		return fdelay;
	}

	void writeSummaryFile(SuiteResult& result)
	{
		std::vector<std::vector<std::string>> rows;
		rows.reserve(result.summary_rows.size());
		for (const auto& [metric, value] : result.summary_rows)
		{
			rows.push_back({metric, value});
		}
		const fs::path path = result.output_dir / "summary.csv";
		writeCsv(path, {"metric", "value"}, rows);
		addGeneratedFile(result, path);
	}

	void writeSourceMappingFile(SuiteResult& result)
	{
		std::vector<std::vector<std::string>> rows;
		rows.reserve(result.source_mappings.size());
		for (const auto& mapping : result.source_mappings)
		{
			rows.push_back({mapping.suite, mapping.upstream_file, mapping.upstream_lines, mapping.local_file, mapping.note});
		}
		const fs::path path = result.output_dir / "source_mapping.csv";
		writeCsv(path, {"suite", "upstream_file", "upstream_lines", "local_file", "note"}, rows);
		addGeneratedFile(result, path);
	}

	void writeReadmeFile(SuiteResult& result)
	{
		const fs::path path = result.output_dir / "README.txt";
		writeTextFile(path, result.readme);
		addGeneratedFile(result, path);
	}

	void writeManifestFile(const SuiteResult& result, const CliOptions& options)
	{
		std::ostringstream stream;
		stream << "{\n";
		stream << "  \"suite\": \"" << jsonEscape(result.suite) << "\",\n";
		stream << "  \"ratio\": " << options.ratio << ",\n";
		stream << "  \"filter_length\": " << options.filter_length << ",\n";
		stream << "  \"seed\": " << options.seed << ",\n";
		stream << "  \"generated_files\": [\n";
		for (size_t i = 0; i < result.generated_files.size(); ++i)
		{
			stream << "    \"" << jsonEscape(result.generated_files[i]) << "\"";
			if (i + 1 < result.generated_files.size())
			{
				stream << ',';
			}
			stream << '\n';
		}
		stream << "  ],\n";
		stream << "  \"source_mappings\": [\n";
		for (size_t i = 0; i < result.source_mappings.size(); ++i)
		{
			const auto& mapping = result.source_mappings[i];
			stream << "    {\n";
			stream << "      \"suite\": \"" << jsonEscape(mapping.suite) << "\",\n";
			stream << "      \"upstream_file\": \"" << jsonEscape(mapping.upstream_file) << "\",\n";
			stream << "      \"upstream_lines\": \"" << jsonEscape(mapping.upstream_lines) << "\",\n";
			stream << "      \"local_file\": \"" << jsonEscape(mapping.local_file) << "\",\n";
			stream << "      \"note\": \"" << jsonEscape(mapping.note) << "\"\n";
			stream << "    }";
			if (i + 1 < result.source_mappings.size())
			{
				stream << ',';
			}
			stream << '\n';
		}
		stream << "  ],\n";
		stream << "  \"deferred_paths\": [\n";
		stream << "    {\n";
		stream << "      \"upstream_file\": \"packages/libfers/src/noise/falpha_branch.cpp\",\n";
		stream << "      \"upstream_lines\": \"154-166\",\n";
		stream << "      \"reason\": \"Documents the DecadeUpsampler caller path but does not execute it in this harness version.\"\n";
		stream << "    }\n";
		stream << "  ]\n";
		stream << "}\n";

		writeTextFile(result.output_dir / "manifest.json", stream.str());
	}

	SuiteResult makeSuiteResult(const std::string& suite, const fs::path& output_dir)
	{
		SuiteResult result;
		result.suite = suite;
		result.output_dir = output_dir;
		ensureDirectory(output_dir);
		return result;
	}

	fers_signal::RadarSignal makeWaveform(const std::string& name, const std::vector<ComplexType>& samples,
										  const RealType sample_rate, const RealType power = 1.0,
										  const RealType carrier = 1.0)
	{
		auto signal = std::make_unique<fers_signal::Signal>();
		signal->load(samples, static_cast<unsigned>(samples.size()), sample_rate);
		return {name, power, carrier, static_cast<RealType>(samples.size()) / sample_rate, std::move(signal)};
	}

	SuiteResult runResamplerSuite(const CliOptions& options)
	{
		ParamGuard guard;
		params::params.reset();
		params::setOversampleRatio(options.ratio);
		params::params.filter_length = options.filter_length;
		params::setRate(1024.0);
		params::setRandomSeed(options.seed);

		SuiteResult result = makeSuiteResult("resampler", options.output_dir);
		result.readme =
			"summary.csv: scalar metrics for the FIR and round-trip checks.\n"
			"fir_coefficients.csv: upstream-equivalent Blackman FIR taps generated from the current ratio and filter length.\n"
			"impulse_response.csv: complex upsampled output for a unit impulse input.\n"
			"step_response.csv: complex upsampled output for a unit step input.\n"
			"tone_roundtrip.csv: tone-by-tone round-trip metrics after upsample+downsample.\n"
			"scenario_roundtrip.csv: multitone, chirp, and low-band noise round-trip metrics.\n"
			"truncation_behavior.csv: output length after downsampling non-divisible oversampled buffers.\n"
			"group_delay_alignment.csv: peak-index alignment checks after round-trip reconstruction.\n";

		result.source_mappings = {
			{"resampler", "packages/libfers/src/signal/dsp_filters.cpp", "71-120", "src/signal/dsp_filters.cpp",
			 "Exact copy of the complex oversample/downsample path under test."},
			{"resampler", "packages/libfers/src/signal/dsp_filters.cpp", "46-65", "src/main.cpp",
			 "Audit-only reimplementation of the upstream Blackman FIR coefficient generator for coefficient export."},
		};

		const auto coeffs = blackmanFirAudit(1.0 / static_cast<RealType>(options.ratio), options.filter_length);
		RealType coeff_sum = std::accumulate(coeffs.begin(), coeffs.end(), 0.0);
		RealType coeff_energy = std::inner_product(coeffs.begin(), coeffs.end(), coeffs.begin(), 0.0);
		const RealType coeff_relative_error =
			options.ratio == 0 ? 0.0 : std::abs(coeff_sum - static_cast<RealType>(options.ratio)) / static_cast<RealType>(options.ratio);
		const RealType estimated_roundtrip_gain =
			options.ratio == 0 ? 0.0 : std::pow(coeff_sum / static_cast<RealType>(options.ratio), 2);

		{
			std::vector<std::vector<std::string>> rows;
			rows.reserve(coeffs.size());
			for (size_t i = 0; i < coeffs.size(); ++i)
			{
				rows.push_back({std::to_string(i), toString(coeffs[i])});
			}
			const fs::path path = result.output_dir / "fir_coefficients.csv";
			writeCsv(path, {"index", "coefficient"}, rows);
			addGeneratedFile(result, path);
		}

		{
			const size_t input_size = 16;
			std::vector<ComplexType> impulse(input_size, ComplexType{0.0, 0.0});
			impulse[0] = ComplexType{1.0, 0.0};
			std::vector<ComplexType> upsampled(input_size * options.ratio);
			fers_signal::upsample(impulse, static_cast<unsigned>(impulse.size()), upsampled);

			std::vector<std::vector<std::string>> rows;
			rows.reserve(upsampled.size());
			for (size_t i = 0; i < upsampled.size(); ++i)
			{
				rows.push_back({std::to_string(i), toString(upsampled[i].real()), toString(upsampled[i].imag()),
								toString(std::abs(upsampled[i]))});
			}
			const fs::path path = result.output_dir / "impulse_response.csv";
			writeCsv(path, {"index", "real", "imag", "magnitude"}, rows);
			addGeneratedFile(result, path);
		}

		{
			const size_t input_size = 16;
			std::vector<ComplexType> step(input_size, ComplexType{1.0, 0.0});
			std::vector<ComplexType> upsampled(input_size * options.ratio);
			fers_signal::upsample(step, static_cast<unsigned>(step.size()), upsampled);

			std::vector<std::vector<std::string>> rows;
			rows.reserve(upsampled.size());
			for (size_t i = 0; i < upsampled.size(); ++i)
			{
				rows.push_back({std::to_string(i), toString(upsampled[i].real()), toString(upsampled[i].imag()),
								toString(std::abs(upsampled[i]))});
			}
			const fs::path path = result.output_dir / "step_response.csv";
			writeCsv(path, {"index", "real", "imag", "magnitude"}, rows);
			addGeneratedFile(result, path);
		}

		{
			const size_t sample_count = 512;
			const std::array<RealType, 5> freqs = {0.01, 0.05, 0.1, 0.18, 0.24};
			std::vector<std::vector<std::string>> rows;
			for (const RealType freq : freqs)
			{
				const auto input = makeTone(sample_count, freq);
				std::vector<ComplexType> upsampled(sample_count * options.ratio);
				fers_signal::upsample(input, static_cast<unsigned>(input.size()), upsampled);
				auto downsampled = fers_signal::downsample(upsampled);
				const auto metrics = compareSignals(input, downsampled, 16);
				rows.push_back({toString(freq, 8), toString(metrics.rms_error), toString(metrics.normalized_correlation),
								toString(metrics.fitted_gain), toString(metrics.fitted_phase)});
			}
			const fs::path path = result.output_dir / "tone_roundtrip.csv";
			writeCsv(path, {"frequency", "rms_error", "normalized_correlation", "fitted_gain", "fitted_phase"}, rows);
			addGeneratedFile(result, path);
		}

		{
			std::mt19937 rng(options.seed);
			std::vector<std::pair<std::string, std::vector<ComplexType>>> scenarios;
			scenarios.emplace_back("multitone", [] {
				std::vector<ComplexType> out(512, ComplexType{0.0, 0.0});
				const auto a = makeTone(512, 0.03);
				const auto b = makeTone(512, 0.07);
				const auto c = makeTone(512, 0.12);
				for (size_t i = 0; i < out.size(); ++i)
				{
					out[i] = (a[i] + b[i] + c[i]) / 3.0;
				}
				return out;
			}());
			scenarios.emplace_back("chirp", makeChirp(512, 0.01, 0.22));
			scenarios.emplace_back("low_band_noise", makeLowBandNoise(512, rng));

			std::vector<std::vector<std::string>> rows;
			for (const auto& [name, input] : scenarios)
			{
				std::vector<ComplexType> upsampled(input.size() * options.ratio);
				fers_signal::upsample(input, static_cast<unsigned>(input.size()), upsampled);
				auto downsampled = fers_signal::downsample(upsampled);
				const auto metrics = compareSignals(input, downsampled, 16);
				rows.push_back({name, toString(metrics.rms_error), toString(metrics.normalized_correlation),
								toString(metrics.fitted_gain), toString(metrics.fitted_phase)});
			}
			const fs::path path = result.output_dir / "scenario_roundtrip.csv";
			writeCsv(path, {"scenario", "rms_error", "normalized_correlation", "fitted_gain", "fitted_phase"}, rows);
			addGeneratedFile(result, path);
		}

		{
			std::vector<std::vector<std::string>> rows;
			for (const unsigned input_size : {63u, 64u, 65u})
			{
				std::vector<ComplexType> input(input_size, ComplexType{0.0, 0.0});
				input[input_size / 2] = {1.0, 0.0};
				std::vector<ComplexType> upsampled(static_cast<size_t>(input_size * options.ratio));
				fers_signal::upsample(input, input_size, upsampled);
				auto downsampled = fers_signal::downsample(upsampled);
				rows.push_back({std::to_string(input_size), std::to_string(upsampled.size()), std::to_string(downsampled.size())});
			}
			const fs::path path = result.output_dir / "truncation_behavior.csv";
			writeCsv(path, {"input_size", "oversampled_size", "downsampled_size"}, rows);
			addGeneratedFile(result, path);
		}

		{
			std::vector<std::vector<std::string>> rows;
			for (const unsigned input_size : {64u, 96u})
			{
				std::vector<ComplexType> input(input_size, ComplexType{0.0, 0.0});
				const unsigned peak_index = input_size / 2;
				input[peak_index] = {1.0, 0.0};
				std::vector<ComplexType> upsampled(static_cast<size_t>(input_size * options.ratio));
				fers_signal::upsample(input, input_size, upsampled);
				auto downsampled = fers_signal::downsample(upsampled);
				const auto max_iter = std::max_element(downsampled.begin(), downsampled.end(),
													   [](const ComplexType& a, const ComplexType& b)
													   { return std::abs(a) < std::abs(b); });
				const unsigned observed_peak = static_cast<unsigned>(std::distance(downsampled.begin(), max_iter));
				rows.push_back({std::to_string(input_size), std::to_string(peak_index), std::to_string(observed_peak),
								toString(static_cast<RealType>(static_cast<int>(observed_peak) - static_cast<int>(peak_index)))});
			}
			const fs::path path = result.output_dir / "group_delay_alignment.csv";
			writeCsv(path, {"input_size", "expected_peak_index", "observed_peak_index", "peak_index_error"}, rows);
			addGeneratedFile(result, path);
		}

		result.summary_rows = {
			{"ratio", toString(options.ratio)},
			{"filter_length", toString(options.filter_length)},
			{"fir_coefficient_count", toString(static_cast<unsigned>(coeffs.size()))},
			{"fir_sum", toString(coeff_sum)},
			{"fir_energy", toString(coeff_energy)},
			{"fir_sum_relative_error", toString(coeff_relative_error)},
			{"estimated_roundtrip_gain", toString(estimated_roundtrip_gain)},
		};

		writeSummaryFile(result);
		writeSourceMappingFile(result);
		writeReadmeFile(result);
		writeManifestFile(result, options);
		return result;
	}

	SuiteResult runRenderSuite(const CliOptions& options)
	{
		ParamGuard guard;
		params::params.reset();
		params::params.filter_length = options.filter_length;
		params::setRate(64.0);
		params::setOversampleRatio(options.ratio);

		SuiteResult result = makeSuiteResult("render", options.output_dir);
		result.readme =
			"summary.csv: scalar metrics for load, interpolation, and oversampled render checks.\n"
			"signal_load.csv: observed versus expected size/rate scaling after Signal::load.\n"
			"constant_render.csv: rendered samples for a constant waveform with fixed amplitude and phase.\n"
			"power_phase_interpolation.csv: midpoint interpolation checks for power and phase.\n"
			"fractional_delay_sweep.csv: fractional-delay wrapping checks against the interpolation filter table.\n"
			"oversampled_render_comparison.csv: ratio-1 render versus oversampled-render-then-downsample metrics.\n"
			"edge_behavior.csv: first and last rendered samples to expose window-edge truncation behavior.\n";

		result.source_mappings = {
			{"render", "packages/libfers/src/signal/radar_signal.cpp", "64-156", "src/signal/radar_signal.cpp",
			 "Exact copy of Signal::load, Signal::render, and RadarSignal::render."},
			{"render", "packages/libfers/src/interpolation/interpolation_filter.cpp", "62-146",
			 "src/interpolation/interpolation_filter.cpp",
			 "Exact copy of the interpolation-filter singleton used during rendering."},
		};

		const unsigned filter_length = params::renderFilterLength();
		const unsigned sample_count = filter_length * 8;
		const std::vector<ComplexType> input(sample_count, ComplexType{1.0, 0.0});
		const auto upsample_coeffs = blackmanFirAudit(1.0 / static_cast<RealType>(options.ratio), filter_length);
		const RealType upsample_coeff_sum = std::accumulate(upsample_coeffs.begin(), upsample_coeffs.end(), 0.0);
		const RealType upsample_dc_gain = upsample_coeff_sum / static_cast<RealType>(options.ratio);

		fers_signal::Signal signal;
		signal.load(input, sample_count, params::rate());

		{
			const fs::path path = result.output_dir / "signal_load.csv";
			writeCsv(path, {"input_samples", "ratio", "observed_rate", "expected_rate", "observed_render_size", "expected_render_size"},
					 {{std::to_string(sample_count), std::to_string(options.ratio), toString(signal.getRate()),
					   toString(params::rate() * options.ratio), std::to_string(sample_count * options.ratio),
					   std::to_string(sample_count * options.ratio)}});
			addGeneratedFile(result, path);
		}

		const auto& interp_filter = interp::InterpFilter::getInstance();

		{
			const RealType power = 4.0;
			const RealType phase = PI / 4.0;
			const std::vector<interp::InterpPoint> points = {{power, 0.0, 0.0, phase}};
			unsigned size = 0;
			const auto data = signal.render(points, size, 0.0);
			const auto filter = interp_filter.getFilter(0.0);
			const auto expected = expectedConstantRender(filter, std::sqrt(power), phase, upsample_dc_gain);

			std::vector<std::vector<std::string>> rows;
			for (size_t i = 0; i < data.size(); ++i)
			{
				rows.push_back({std::to_string(i), toString(data[i].real()), toString(data[i].imag()), toString(std::abs(data[i]))});
			}
			const fs::path path = result.output_dir / "constant_render.csv";
			writeCsv(path, {"index", "real", "imag", "magnitude"}, rows);
			addGeneratedFile(result, path);

			result.summary_rows.push_back({"constant_expected_real", toString(expected.real())});
			result.summary_rows.push_back({"constant_expected_imag", toString(expected.imag())});
			result.summary_rows.push_back({"constant_observed_center_real", toString(data[filter_length].real())});
			result.summary_rows.push_back({"constant_observed_center_imag", toString(data[filter_length].imag())});
			result.summary_rows.push_back({"upsample_dc_gain", toString(upsample_dc_gain)});
		}

		{
			const RealType power_a = 1.0;
			const RealType power_b = 9.0;
			const RealType phase_a = 0.0;
			const RealType phase_b = PI / 2.0;
			const std::vector<interp::InterpPoint> points = {{power_a, 0.0, 0.0, phase_a},
															 {power_b, 2.0 * filter_length, 0.0, phase_b}};

			unsigned size = 0;
			const auto data = signal.render(points, size, 0.0);
			const auto filter = interp_filter.getFilter(0.0);

			std::vector<std::vector<std::string>> rows;
			for (unsigned index : {filter_length / 2, filter_length, filter_length + filter_length / 2, filter_length * 2})
			{
				if (index >= data.size())
				{
					continue;
				}

				const RealType sample_time = static_cast<RealType>(index) / signal.getRate();
				const RealType bw = sample_time / (2.0 * filter_length / signal.getRate());
				const RealType amplitude = std::lerp(std::sqrt(power_a), std::sqrt(power_b), std::clamp(bw, 0.0, 1.0));
				const RealType phase = std::lerp(phase_a, phase_b, std::clamp(bw, 0.0, 1.0));
				const auto expected = expectedConstantRender(filter, amplitude, phase, upsample_dc_gain);
				rows.push_back({std::to_string(index), toString(sample_time), toString(expected.real()), toString(expected.imag()),
								toString(data[index].real()), toString(data[index].imag())});
			}
			const fs::path path = result.output_dir / "power_phase_interpolation.csv";
			writeCsv(path, {"index", "sample_time", "expected_real", "expected_imag", "observed_real", "observed_imag"}, rows);
			addGeneratedFile(result, path);
		}

		{
			const std::vector<interp::InterpPoint> points = {{1.0, 0.0, 0.0, 0.0}};
			std::vector<std::vector<std::string>> rows;
			for (const RealType frac_delay : {-0.95, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 0.95})
			{
				unsigned size = 0;
				const auto data = signal.render(points, size, frac_delay);
				const RealType effective_delay = wrappedDelayFromFracWinDelay(frac_delay);
				const auto filter = interp_filter.getFilter(effective_delay);
				const auto expected = expectedConstantRender(filter, 1.0, 0.0, upsample_dc_gain);
				rows.push_back({toString(frac_delay, 8), toString(effective_delay, 8), toString(expected.real()),
								toString(data[filter_length].real()), toString(data[filter_length].imag())});
			}
			const fs::path path = result.output_dir / "fractional_delay_sweep.csv";
			writeCsv(path, {"frac_win_delay", "effective_delay", "expected_real", "observed_real", "observed_imag"}, rows);
			addGeneratedFile(result, path);
		}

		{
			params::setOversampleRatio(1);
			params::setRate(64.0);
			const auto reference_waveform = makeTone(256, 0.07);
			fers_signal::Signal reference_signal;
			reference_signal.load(reference_waveform, static_cast<unsigned>(reference_waveform.size()), params::rate());
			const std::vector<interp::InterpPoint> points = {{1.0, 0.0, 0.0, 0.0}};
			unsigned reference_size = 0;
			const auto reference_render = reference_signal.render(points, reference_size, 0.15);

			params::setOversampleRatio(options.ratio);
			fers_signal::Signal oversampled_signal;
			oversampled_signal.load(reference_waveform, static_cast<unsigned>(reference_waveform.size()), 64.0);
			unsigned oversampled_size = 0;
			auto oversampled_render = oversampled_signal.render(points, oversampled_size, 0.15);
			auto downsampled_render = fers_signal::downsample(oversampled_render);
			const auto metrics = compareSignals(reference_render, downsampled_render, 16);

			const fs::path path = result.output_dir / "oversampled_render_comparison.csv";
			writeCsv(path, {"scenario", "rms_error", "normalized_correlation", "fitted_gain", "fitted_phase"},
					 {{"tone_constant_delay", toString(metrics.rms_error), toString(metrics.normalized_correlation),
					   toString(metrics.fitted_gain), toString(metrics.fitted_phase)}});
			addGeneratedFile(result, path);

			result.summary_rows.push_back({"oversampled_render_rms_error", toString(metrics.rms_error)});
			result.summary_rows.push_back({"oversampled_render_correlation", toString(metrics.normalized_correlation)});
		}

		{
			params::setOversampleRatio(options.ratio);
			fers_signal::Signal edge_signal;
			edge_signal.load(input, sample_count, 64.0);
			const std::vector<interp::InterpPoint> points = {{1.0, 0.0, 0.0, 0.0}};
			unsigned size = 0;
			const auto data = edge_signal.render(points, size, 0.0);

			std::vector<std::vector<std::string>> rows;
			for (size_t i = 0; i < std::min<size_t>(16, data.size()); ++i)
			{
				rows.push_back({"head", std::to_string(i), toString(data[i].real()), toString(data[i].imag()),
								toString(std::abs(data[i]))});
			}
			for (size_t i = data.size() > 16 ? data.size() - 16 : 0; i < data.size(); ++i)
			{
				rows.push_back({"tail", std::to_string(i), toString(data[i].real()), toString(data[i].imag()),
								toString(std::abs(data[i]))});
			}
			const fs::path path = result.output_dir / "edge_behavior.csv";
			writeCsv(path, {"region", "index", "real", "imag", "magnitude"}, rows);
			addGeneratedFile(result, path);
		}

		result.summary_rows.insert(result.summary_rows.begin(),
								   {{"ratio", toString(options.ratio)}, {"filter_length", toString(options.filter_length)}});
		writeSummaryFile(result);
		writeSourceMappingFile(result);
		writeReadmeFile(result);
		writeManifestFile(result, options);
		return result;
	}

	SuiteResult runPipelineSuite(const CliOptions& options)
	{
		ParamGuard guard;
		params::params.reset();
		params::setOversampleRatio(options.ratio);
		params::params.filter_length = options.filter_length;
		params::setRate(64.0);
		params::setRandomSeed(options.seed);

		SuiteResult result = makeSuiteResult("pipeline", options.output_dir);
		result.readme =
			"summary.csv: scalar metrics for renderWindow, thermal noise, and final decimation.\n"
			"response_overlap.csv: window overlap and clipping metrics for multiple responses.\n"
			"window_trace.csv: rendered receive-window samples after response accumulation.\n"
			"thermal_noise_stats.csv: empirical means and variances against the implemented bandwidth formula.\n"
			"downsampling_quantization.csv: final downsample-plus-normalize metrics.\n"
			"adc_histogram_bits_*.csv: quantized code histograms for representative ADC depths.\n";

		result.source_mappings = {
			{"pipeline", "packages/libfers/src/serial/response.cpp", "27-33", "src/serial/response.cpp",
			 "Exact copy of renderBinary delegation."},
			{"pipeline", "packages/libfers/src/processing/signal_processor.cpp", "59-154", "src/processing/signal_processor.cpp",
			 "Source-faithful renderWindow, thermal noise, and quantization behavior."},
			{"pipeline", "packages/libfers/src/processing/finalizer_pipeline.cpp", "113-120",
			 "src/processing/finalizer_pipeline.cpp", "Exact copy of the oversampling-aware final decimation step."},
		};

		const auto waveform_samples = makeTone(16, 0.05);
		auto wave = makeWaveform("audit_wave", waveform_samples, 64.0);

		std::vector<std::unique_ptr<serial::Response>> responses;
		{
			auto before = std::make_unique<serial::Response>(&wave);
			before->addInterpPoint({1.0, 0.875, 0.0, 0.0});
			before->addInterpPoint({1.0, 1.125, 0.0, 0.0});
			responses.push_back(std::move(before));

			auto inside = std::make_unique<serial::Response>(&wave);
			inside->addInterpPoint({0.5, 1.125, 0.0, PI / 8.0});
			inside->addInterpPoint({0.5, 1.375, 0.0, PI / 8.0});
			responses.push_back(std::move(inside));

			auto after = std::make_unique<serial::Response>(&wave);
			after->addInterpPoint({0.75, 1.375, 0.0, PI / 4.0});
			after->addInterpPoint({0.75, 1.625, 0.0, PI / 4.0});
			responses.push_back(std::move(after));
		}

		const RealType window_start = 1.0;
		const RealType window_length = 0.5;
		const unsigned window_size = static_cast<unsigned>(std::ceil(window_length * params::rate() * params::oversampleRatio()));
		std::vector<ComplexType> window(window_size, ComplexType{0.0, 0.0});
		processing::renderWindow(window, window_length, window_start, 0.0, responses);

		{
			std::vector<std::vector<std::string>> overlap_rows;
			for (const auto& response : responses)
			{
				RealType rate = 0.0;
				unsigned size = 0;
				const auto rendered = response->renderBinary(rate, size, 0.0);
				const int start_sample = static_cast<int>(std::round(params::rate() * params::oversampleRatio() *
																 (response->startTime() - window_start)));
				const unsigned clipped_prefix = start_sample < 0 ? static_cast<unsigned>(-start_sample) : 0u;
				const unsigned clipped_suffix = start_sample + static_cast<int>(size) > static_cast<int>(window_size) ?
					static_cast<unsigned>(start_sample + static_cast<int>(size) - static_cast<int>(window_size)) :
					0u;
				RealType energy = 0.0;
				for (const auto& sample : rendered)
				{
					energy += std::norm(sample);
				}

				overlap_rows.push_back({toString(response->startTime()), toString(response->endTime()), std::to_string(size),
										std::to_string(clipped_prefix), std::to_string(clipped_suffix), toString(energy)});
			}
			const fs::path path = result.output_dir / "response_overlap.csv";
			writeCsv(path, {"start_time", "end_time", "render_size", "clipped_prefix_samples", "clipped_suffix_samples",
							"render_energy"},
					 overlap_rows);
			addGeneratedFile(result, path);
		}

		{
			std::vector<std::vector<std::string>> trace_rows;
			trace_rows.reserve(window.size());
			for (size_t i = 0; i < window.size(); ++i)
			{
				trace_rows.push_back({std::to_string(i), toString(window[i].real()), toString(window[i].imag()),
									  toString(std::abs(window[i]))});
			}
			const fs::path path = result.output_dir / "window_trace.csv";
			writeCsv(path, {"index", "real", "imag", "magnitude"}, trace_rows);
			addGeneratedFile(result, path);
		}

		{
			std::vector<std::vector<std::string>> rows;
			for (const RealType temperature : {0.0, 50.0, 300.0})
			{
				std::vector<ComplexType> noise_buffer(50000, ComplexType{0.0, 0.0});
				std::mt19937 rng(options.seed + static_cast<unsigned>(temperature));
				processing::applyThermalNoise(noise_buffer, temperature, rng);

				RealType mean_real = 0.0;
				RealType mean_imag = 0.0;
				for (const auto& sample : noise_buffer)
				{
					mean_real += sample.real();
					mean_imag += sample.imag();
				}
				mean_real /= static_cast<RealType>(noise_buffer.size());
				mean_imag /= static_cast<RealType>(noise_buffer.size());

				RealType var_real = 0.0;
				RealType var_imag = 0.0;
				for (const auto& sample : noise_buffer)
				{
					var_real += std::pow(sample.real() - mean_real, 2);
					var_imag += std::pow(sample.imag() - mean_imag, 2);
				}
				var_real /= static_cast<RealType>(noise_buffer.size());
				var_imag /= static_cast<RealType>(noise_buffer.size());

				const RealType bandwidth = params::rate() / (2.0 * params::oversampleRatio());
				const RealType expected_power = params::boltzmannK() * temperature * bandwidth / 2.0;
				rows.push_back({toString(temperature), toString(mean_real), toString(mean_imag), toString(var_real),
								toString(var_imag), toString(expected_power)});
			}
			const fs::path path = result.output_dir / "thermal_noise_stats.csv";
			writeCsv(path, {"temperature", "mean_real", "mean_imag", "variance_real", "variance_imag",
							"expected_per_channel_power"},
					 rows);
			addGeneratedFile(result, path);
		}

		{
			params::setAdcBits(0);
			const size_t sample_count = 64;
			const auto baseband = makeTone(sample_count, 0.05);
			std::vector<ComplexType> oversampled(sample_count * params::oversampleRatio());
			fers_signal::upsample(baseband, static_cast<unsigned>(baseband.size()), oversampled);
			const RealType fullscale = processing::pipeline::applyDownsamplingAndQuantization(oversampled);
			const auto metrics = compareSignals(baseband, oversampled, 10);

			const fs::path path = result.output_dir / "downsampling_quantization.csv";
			writeCsv(path, {"scenario", "final_size", "fullscale", "rms_error", "normalized_correlation"},
					 {{"normalize_only", std::to_string(oversampled.size()), toString(fullscale), toString(metrics.rms_error),
					   toString(metrics.normalized_correlation)}});
			addGeneratedFile(result, path);

			result.summary_rows.push_back({"downsampled_fullscale", toString(fullscale)});
			result.summary_rows.push_back({"downsampled_correlation", toString(metrics.normalized_correlation)});
		}

		for (const unsigned bits : {3u, 8u})
		{
			params::setAdcBits(bits);
			std::vector<ComplexType> adc_window;
			adc_window.reserve(256);
			for (int i = 0; i < 256; ++i)
			{
				const RealType x = -1.0 + 2.0 * static_cast<RealType>(i) / 255.0;
				adc_window.emplace_back(x, 0.5 * x);
			}

			const RealType fullscale = processing::quantizeAndScaleWindow(adc_window);
			(void)fullscale;

			std::map<std::string, unsigned> histogram;
			for (const auto& sample : adc_window)
			{
				++histogram[toString(sample.real(), 6)];
			}

			std::vector<std::vector<std::string>> rows;
			for (const auto& [level, count] : histogram)
			{
				rows.push_back({level, std::to_string(count)});
			}
			const fs::path path = result.output_dir / ("adc_histogram_bits_" + std::to_string(bits) + ".csv");
			writeCsv(path, {"quantized_real_level", "count"}, rows);
			addGeneratedFile(result, path);
		}

		result.summary_rows.insert(result.summary_rows.begin(),
								   {{"ratio", toString(options.ratio)}, {"filter_length", toString(options.filter_length)}});
		writeSummaryFile(result);
		writeSourceMappingFile(result);
		writeReadmeFile(result);
		writeManifestFile(result, options);
		return result;
	}

	RealType quantizedPrf(const RealType rate, const unsigned ratio, const RealType requested_prf)
	{
		const RealType oversampled_rate = rate * ratio;
		return 1.0 / (std::floor(oversampled_rate / requested_prf) / oversampled_rate);
	}

	RealType quantizedSkip(const RealType rate, const unsigned ratio, const RealType skip)
	{
		const RealType oversampled_rate = rate * ratio;
		return std::floor(oversampled_rate * skip) / oversampled_rate;
	}

	size_t cwBufferSize(const RealType start_time, const RealType end_time, const RealType rate, const unsigned ratio)
	{
		const RealType dt_sim = 1.0 / (rate * ratio);
		return static_cast<size_t>(std::ceil((end_time - start_time) / dt_sim));
	}

	SuiteResult runClockingSuite(const CliOptions& options)
	{
		SuiteResult result = makeSuiteResult("clocking", options.output_dir);
		result.readme =
			"summary.csv: scalar maxima for the caller-formula probes.\n"
			"transmitter_prf_quantization.csv: requested versus quantized transmitter PRFs.\n"
			"receiver_window_quantization.csv: quantized receiver window PRFs and skips.\n"
			"cw_timestep_and_buffer.csv: dt_sim and CW buffer sizes for representative durations.\n"
			"source_mapping.csv: caller-to-formula mapping for non-sample-transforming oversampling users.\n";

		result.source_mappings = {
			{"clocking", "packages/libfers/src/radar/transmitter.cpp", "19-23", "src/main.cpp",
			 "Formula probe mirroring transmitter PRF quantization."},
			{"clocking", "packages/libfers/src/radar/receiver.cpp", "90-96", "src/main.cpp",
			 "Formula probe mirroring receiver window PRF and skip quantization."},
			{"clocking", "packages/libfers/src/core/sim_threading.cpp", "103-110", "src/main.cpp",
			 "Formula probe mirroring the oversampled CW timestep."},
			{"clocking", "packages/libfers/src/serial/xml_parser_utils.cpp", "1071-1081", "src/main.cpp",
			 "Formula probe mirroring CW buffer preallocation from parsed parameters."},
			{"clocking", "packages/libfers/src/serial/json_serializer.cpp", "1583-1594", "src/main.cpp",
			 "Formula probe mirroring CW buffer preallocation from hydrated JSON state."},
		};

		RealType max_prf_error = 0.0;
		{
			std::vector<std::vector<std::string>> rows;
			for (const RealType rate : {1.0e3, 1.0e6})
			{
				for (const unsigned ratio : {1u, 2u, 4u, options.ratio})
				{
					for (const RealType requested_prf : {7.0, 33.0, 123.45, 499.9})
					{
						const RealType effective = quantizedPrf(rate, ratio, requested_prf);
						const RealType error = effective - requested_prf;
						max_prf_error = std::max(max_prf_error, std::abs(error));
						rows.push_back({toString(rate), std::to_string(ratio), toString(rate * ratio), toString(requested_prf),
										toString(effective), toString(error)});
					}
				}
			}
			const fs::path path = result.output_dir / "transmitter_prf_quantization.csv";
			writeCsv(path, {"rate", "ratio", "oversampled_rate", "requested_prf", "effective_prf", "error"}, rows);
			addGeneratedFile(result, path);
		}

		RealType max_skip_error = 0.0;
		{
			std::vector<std::vector<std::string>> rows;
			for (const RealType rate : {1.0e3, 1.0e6})
			{
				for (const unsigned ratio : {1u, 2u, 4u, options.ratio})
				{
					for (const RealType requested_prf : {7.0, 33.0, 123.45})
					{
						for (const RealType skip : {0.0, 1.0e-6, 0.00037, 0.0013})
						{
							const RealType effective_prf = quantizedPrf(rate, ratio, requested_prf);
							const RealType effective_skip = quantizedSkip(rate, ratio, skip);
							max_skip_error = std::max(max_skip_error, std::abs(effective_skip - skip));
							rows.push_back({toString(rate), std::to_string(ratio), toString(requested_prf), toString(effective_prf),
											toString(skip), toString(effective_skip), toString(effective_skip - skip)});
						}
					}
				}
			}
			const fs::path path = result.output_dir / "receiver_window_quantization.csv";
			writeCsv(path, {"rate", "ratio", "requested_prf", "effective_prf", "requested_skip", "effective_skip", "skip_error"},
					 rows);
			addGeneratedFile(result, path);
		}

		{
			std::vector<std::vector<std::string>> rows;
			for (const RealType rate : {1.0e3, 1.0e6})
			{
				for (const unsigned ratio : {1u, 2u, 4u, options.ratio})
				{
					for (const auto [start_time, end_time] :
						 {std::pair{0.0, 1.0}, std::pair{0.0, 1.001}, std::pair{0.125, 0.125123}})
					{
						const RealType dt_sim = 1.0 / (rate * ratio);
						rows.push_back({toString(rate), std::to_string(ratio), toString(start_time), toString(end_time),
										toString(dt_sim), std::to_string(cwBufferSize(start_time, end_time, rate, ratio))});
					}
				}
			}
			const fs::path path = result.output_dir / "cw_timestep_and_buffer.csv";
			writeCsv(path, {"rate", "ratio", "start_time", "end_time", "dt_sim", "cw_buffer_size"}, rows);
			addGeneratedFile(result, path);
		}

		result.summary_rows = {{"max_abs_prf_error", toString(max_prf_error)}, {"max_abs_skip_error", toString(max_skip_error)}};
		writeSummaryFile(result);
		writeSourceMappingFile(result);
		writeReadmeFile(result);
		writeManifestFile(result, options);
		return result;
	}

	SuiteResult runSuite(const CliOptions& options)
	{
		if (options.suite == "resampler")
		{
			return runResamplerSuite(options);
		}
		if (options.suite == "render")
		{
			return runRenderSuite(options);
		}
		if (options.suite == "pipeline")
		{
			return runPipelineSuite(options);
		}
		if (options.suite == "clocking")
		{
			return runClockingSuite(options);
		}
		throw std::runtime_error("Unknown suite: " + options.suite);
	}

	int runChildSuite(const CliOptions& parent_options, const std::string& suite, const unsigned ratio,
					  const unsigned filter_length, const fs::path& output_dir)
	{
		std::vector<std::string> args = {"suite",
										 suite,
										 "--ratio",
										 std::to_string(ratio),
										 "--filter-length",
										 std::to_string(filter_length),
										 "--seed",
										 std::to_string(parent_options.seed),
										 "--output-dir",
										 output_dir.string()};

		std::ostringstream command;
		command << shellQuote(parent_options.executable_path.string());
		for (const auto& arg : args)
		{
			command << ' ' << shellQuote(arg);
		}
		return std::system(command.str().c_str());
	}

	int runAll(const CliOptions& options)
	{
		const std::array<std::string_view, 4> suites = {"resampler", "render", "pipeline", "clocking"};
		for (const auto suite : suites)
		{
			const fs::path child_output =
				options.output_dir / ("ratio_" + std::to_string(options.ratio)) / ("filter_" + std::to_string(options.filter_length)) /
				std::string(suite);
			const int status = runChildSuite(options, std::string(suite), options.ratio, options.filter_length, child_output);
			if (status != 0)
			{
				return status;
			}
		}
		return 0;
	}

	int runMatrix(const CliOptions& options)
	{
		for (const unsigned ratio : options.ratios)
		{
			for (const unsigned filter_length : options.filter_lengths)
			{
				const std::array<std::string_view, 4> suites = {"resampler", "render", "pipeline", "clocking"};
				for (const auto suite : suites)
				{
					const fs::path child_output =
						options.output_dir / ("ratio_" + std::to_string(ratio)) / ("filter_" + std::to_string(filter_length)) /
						std::string(suite);
					const int status = runChildSuite(options, std::string(suite), ratio, filter_length, child_output);
					if (status != 0)
					{
						return status;
					}
				}
			}
		}
		return 0;
	}

	std::string usage()
	{
		return "Usage:\n"
			   "  oversampling_audit all [--ratio N] [--filter-length N] [--seed N] [--output-dir PATH]\n"
			   "  oversampling_audit suite <resampler|render|pipeline|clocking> [--ratio N] [--filter-length N] "
			   "[--seed N] [--output-dir PATH]\n"
			   "  oversampling_audit matrix [--ratios csv] [--filter-lengths csv] [--seed N] [--output-dir PATH]\n";
	}
}

int main(const int argc, char* argv[])
{
	try
	{
		logging::logger.setLevel(logging::Level::WARNING);
		const CliOptions options = parseCli(argc, argv);

		switch (options.command)
		{
		case CommandType::Suite:
			runSuite(options);
			return 0;
		case CommandType::All:
			return runAll(options);
		case CommandType::Matrix:
			return runMatrix(options);
		}
	}
	catch (const std::exception& error)
	{
		std::cerr << error.what() << '\n' << usage();
		return 1;
	}

	return 1;
}
