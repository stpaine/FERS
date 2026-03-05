cmake_minimum_required(VERSION 3.22)

project(FERS
	VERSION 1.0.0
	DESCRIPTION "The Flexible Extensible Radar Simulator"
	LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# --- Project-wide Settings ---
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON) # For simpler DLL exports

# --- Options and Build Type ---
include(FersOptions)

# --- Compiler Warnings and Flags ---
include(FersCompilerWarnings)

# --- Testing ---
if (FERS_BUILD_TESTS)
	enable_testing()
endif()

if(FERS_ENABLE_COVERAGE)
	if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
		message(STATUS "Enabling Code Coverage")
		add_compile_options(-fprofile-arcs -ftest-coverage)
		add_link_options(--coverage)
	else()
		message(WARNING "Code coverage is only supported for GCC and Clang")
	endif()
endif()

# --- macOS Specific Configuration ---
if (APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	message(STATUS "macOS with a modern toolchain (e.g., Homebrew LLVM) is recommended.")
	message(STATUS "If you are not using the system compiler, please configure with:")
	message(STATUS "  CC=/path/to/clang CXX=/path/to/clang++ cmake ..")
	# link_directories(/usr/local/opt/llvm/lib/c++)
endif ()

# Add the subdirectory containing the C++ packages.
add_subdirectory(packages)

# --- Documentation ---
add_subdirectory(docs)
