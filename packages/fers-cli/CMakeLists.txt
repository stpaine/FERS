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

# Link against libfers. CMake handles transitive dependencies and include paths.
target_link_libraries(fers-cli PRIVATE libfers)

# Install the executable
install(TARGETS fers-cli DESTINATION bin)
