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
			   signal/test_dsp_filters.cpp
			   signal/test_radar_signal.cpp
)

target_include_directories(libfers_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)

target_link_libraries(libfers_tests PRIVATE Catch2::Catch2WithMain fers::fers)

if(COMMAND apply_fers_warnings)
	apply_fers_warnings(libfers_tests)
endif()

catch_discover_tests(libfers_tests)
