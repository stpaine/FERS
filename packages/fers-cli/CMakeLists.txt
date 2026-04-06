project(fers-cli)

# Create the internal library target that holds the unit-testable CLI logic.
add_library(fers-cli-core STATIC)

apply_fers_warnings(fers-cli-core)

target_sources(fers-cli-core PRIVATE
			   src/arg_parser.cpp
			   src/cli_paths.cpp
)

target_include_directories(fers-cli-core PUBLIC
						   ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(fers-cli-core PUBLIC fers::fers)

# Create the executable target
add_executable(fers-cli)

# Apply project-wide compiler warnings and flags
apply_fers_warnings(fers-cli)

target_sources(fers-cli PRIVATE
			   src/main.cpp
)

target_link_libraries(fers-cli PRIVATE fers-cli-core)

# --- Cross-Platform Relocatable RPATH ---
set_target_properties(fers-cli PROPERTIES
					  BUILD_WITH_INSTALL_RPATH FALSE
					  INSTALL_RPATH_USE_LINK_PATH TRUE
)

if (APPLE)
	set_target_properties(fers-cli PROPERTIES INSTALL_RPATH "@executable_path/../lib")
elseif (UNIX)
	set_target_properties(fers-cli PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
endif ()

# Install the executable
install(TARGETS fers-cli DESTINATION bin)

if (FERS_BUILD_TESTS)
	add_subdirectory(tests)
endif ()
