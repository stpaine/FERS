project(fers-cli)

# Create the executable target
add_executable(fers-cli)

# Apply project-wide compiler warnings and flags
apply_fers_warnings(fers-cli)

# Add its own source files
target_sources(fers-cli PRIVATE
			   src/main.cpp
			   src/arg_parser.cpp
)

# Add include directory for its own headers (arg_parser.h)
target_include_directories(fers-cli PRIVATE
						   ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Link against libfers.
target_link_libraries(fers-cli PRIVATE fers::fers_static)

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
