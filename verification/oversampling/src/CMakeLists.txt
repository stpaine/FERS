add_executable(oversampling_audit
	main.cpp
	core/logging.cpp
	interpolation/interpolation_filter.cpp
	processing/finalizer_pipeline.cpp
	processing/signal_processor.cpp
	serial/response.cpp
	signal/dsp_filters.cpp
	signal/radar_signal.cpp
)

target_include_directories(oversampling_audit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

