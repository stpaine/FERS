include(FetchContent)

FetchContent_Declare(
		Catch2
		GIT_REPOSITORY https://github.com/catchorg/Catch2.git
		GIT_TAG        v3.13.0
)
FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)

add_executable(libfers_tests
			   math/test_coord.cpp
			   math/test_geometry_ops.cpp
			   math/test_path.cpp
			   math/test_rotation_path.cpp
			   interpolation/test_interpolation_set.cpp
			   interpolation/test_interpolation_filter.cpp
			   noise/test_noise_generators.cpp
			   noise/test_falpha_branch.cpp
			   antenna/test_antenna_factory.cpp
			   core/test_config.cpp
			   core/test_fers_context.cpp
			   core/test_logging.cpp
			   core/test_parameters.cpp
			   core/test_portable_utils.cpp
			   core/test_sim_events.cpp
			   core/test_sim_id.cpp
			   core/test_thread_pool.cpp
			   signal/test_dsp_filters.cpp
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
)

target_include_directories(libfers_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)

target_link_libraries(libfers_tests PRIVATE Catch2::Catch2WithMain fers::fers)

if(COMMAND apply_fers_warnings)
	apply_fers_warnings(libfers_tests)
endif()

catch_discover_tests(libfers_tests)
