find_package(Catch2 CONFIG REQUIRED)

include(CTest)
include(Catch)

add_executable(fers_cli_tests
			   test_arg_parser.cpp
			   test_cli_paths.cpp
)

target_link_libraries(fers_cli_tests PRIVATE
					  Catch2::Catch2WithMain
					  fers-cli-core
)

if (COMMAND apply_fers_warnings)
	apply_fers_warnings(fers_cli_tests)
endif ()

catch_discover_tests(fers_cli_tests
					 DISCOVERY_MODE PRE_TEST
)
