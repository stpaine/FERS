find_package(Catch2 CONFIG REQUIRED)

include(CTest)
include(Catch)

add_executable(libfers_tests
               test_main.cpp
			   math/test_coord.cpp
			   math/test_geometry_ops.cpp
			   math/test_path.cpp
			   math/test_rotation_path.cpp
			   interpolation/test_interpolation_set.cpp
			   interpolation/test_interpolation_filter.cpp
			   noise/test_noise_generators.cpp
			   noise/test_falpha_branch.cpp
			   antenna/test_antenna_factory.cpp
			   antenna/test_antenna_factory_io.cpp
			   api/test_api_contract.cpp
			   api/test_api_scenario_io.cpp
			   api/test_api_runtime.cpp
			   api/test_api_paths.cpp
			   api/test_api_preview.cpp
			   core/test_config.cpp
			   core/test_fers_context.cpp
			   core/test_logging.cpp
			   core/test_memory_projection.cpp
			   core/test_parameters.cpp
			   core/test_sim_threading.cpp
			   core/test_portable_utils.cpp
			   core/test_sim_events.cpp
			   core/test_sim_id.cpp
			   core/test_thread_pool.cpp
			   core/test_world.cpp
			   signal/test_dsp_filters.cpp
			   signal/test_if_resampler.cpp
			   signal/test_radar_signal.cpp
			   timing/test_prototype_timing.cpp
			   timing/test_timing.cpp
			   radar/test_object.cpp
			   radar/test_platform.cpp
			   radar/test_schedule_period.cpp
			   radar/test_radar_obj.cpp
			   radar/test_transmitter.cpp
			   radar/test_target.cpp
			   radar/test_receiver.cpp
			   processing/test_signal_processor.cpp
			   processing/test_finalizer_pipeline_timing.cpp
			   processing/test_finalizer_pipeline_interference.cpp
			   processing/test_finalizer_pipeline_io.cpp
			   processing/test_finalizer.cpp
			   serial/test_response.cpp
			   serial/test_libxml_wrapper.cpp
			   serial/test_hdf5_handler.cpp
			   serial/test_waveform_factory.cpp
			   serial/test_xml_parser.cpp
			   serial/test_xml_parser_utils.cpp
			   serial/test_rotation_warning_utils.cpp
			   serial/test_json_serializer_stress.cpp
			   serial/test_json_serializer.cpp
			   serial/test_xml_serializer.cpp
			   serial/test_kml_generator.cpp
			   serial/test_kml_generator_utils.cpp
			   simulation/test_channel_model_helpers.cpp
			   simulation/test_channel_model_direct.cpp
			   simulation/test_channel_model_reflected.cpp
			   simulation/test_channel_model_response.cpp
			   simulation/test_channel_model_preview.cpp
)

target_include_directories(libfers_tests PRIVATE
						   ${CMAKE_CURRENT_SOURCE_DIR}/../src
						   ${HDF5_INCLUDE_DIRS})

# --- macOS HDF5 ODR Violation ---
# When HDF5 is built as a static library (vcpkg default), linking it into both
# a shared libfers and the test executable creates two separate HDF5 global states.
# By forcing the tests to link against the static libfers, we ensure HDF5 is
# statically linked exactly once into the final test binary.
if (TARGET fers::fers_static)
	set(FERS_TEST_TARGET fers::fers_static)
else ()
	set(FERS_TEST_TARGET fers::fers)
endif ()

target_link_libraries(libfers_tests PRIVATE
					  Catch2::Catch2WithMain
					  ${FERS_TEST_TARGET}
					  LibXml2::LibXml2
					  # Using this format instead of hdf5::hdf5/hdf5::hdf5_hl due to macOS linking issues
					  ${HDF5_C_LIBRARIES}
					  ${HDF5_HL_LIBRARIES}
)

if (COMMAND apply_fers_warnings)
	apply_fers_warnings(libfers_tests)
endif ()

target_compile_definitions(libfers_tests PRIVATE
						   FERS_TEST_EXPECTED_VERSION="${CMAKE_PROJECT_VERSION}")

catch_discover_tests(libfers_tests
					 DISCOVERY_MODE PRE_TEST
)
