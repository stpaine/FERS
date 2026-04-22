cmake_minimum_required(VERSION 3.22)

set(FERS_VERSION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/version.txt")
if (NOT EXISTS "${FERS_VERSION_FILE}")
	message(FATAL_ERROR "Missing version file: ${FERS_VERSION_FILE}")
endif ()

file(READ "${FERS_VERSION_FILE}" FERS_PROJECT_VERSION)
string(STRIP "${FERS_PROJECT_VERSION}" FERS_PROJECT_VERSION)

if (NOT FERS_PROJECT_VERSION MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$")
	message(FATAL_ERROR "Invalid FERS version '${FERS_PROJECT_VERSION}' in ${FERS_VERSION_FILE}")
endif ()

# --- Enforce macOS Deployment Target ---
if (APPLE)
	if (NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET)
		set(CMAKE_OSX_DEPLOYMENT_TARGET "14.0" CACHE STRING "Minimum macOS version")
	endif ()
	# Ensure vcpkg also respects this deployment target
	set(VCPKG_OSX_DEPLOYMENT_TARGET ${CMAKE_OSX_DEPLOYMENT_TARGET} CACHE STRING "")
endif ()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# --- Options and Build Type ---
include(FersOptions)

if (NOT FERS_DOCS_ONLY)
	# Auto-detect vcpkg toolchain before project() so dependency resolution works in CLion and CLI.
	# Note: when building via `cargo build` (the Tauri path), CMAKE_TOOLCHAIN_FILE is always passed
	# explicitly by build.rs, so this block is effectively dead in that path. It is retained here
	# for developers who invoke CMake directly (e.g. cmake --preset=release from the repo root).
	if (NOT DEFINED VCPKG_TARGET_TRIPLET)
		if (DEFINED ENV{VCPKG_TARGET_TRIPLET})
			set(VCPKG_TARGET_TRIPLET "$ENV{VCPKG_TARGET_TRIPLET}" CACHE STRING "vcpkg target triplet")
		elseif (CMAKE_HOST_WIN32)
			if (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "^(ARM64|AARCH64|aarch64)")
				set(VCPKG_TARGET_TRIPLET "arm64-windows-static-md" CACHE STRING "vcpkg target triplet")
			else ()
				set(VCPKG_TARGET_TRIPLET "x64-windows-static-md" CACHE STRING "vcpkg target triplet")
			endif ()
		endif ()
	endif ()

	if (NOT DEFINED CMAKE_TOOLCHAIN_FILE)
		set(_fers_vcpkg_candidates)

		if (DEFINED ENV{VCPKG_ROOT})
			list(APPEND _fers_vcpkg_candidates "$ENV{VCPKG_ROOT}")
		endif ()

		if (DEFINED ENV{HOME})
			list(APPEND _fers_vcpkg_candidates "$ENV{HOME}/vcpkg")
		endif ()

		if (DEFINED ENV{USERPROFILE})
			list(APPEND _fers_vcpkg_candidates "$ENV{USERPROFILE}/vcpkg")
		endif ()

		list(APPEND _fers_vcpkg_candidates "${CMAKE_SOURCE_DIR}/vcpkg")

		find_program(_fers_vcpkg_exe NAMES vcpkg)
		if (_fers_vcpkg_exe)
			cmake_path(GET _fers_vcpkg_exe PARENT_PATH _fers_vcpkg_bin_dir)
			cmake_path(GET _fers_vcpkg_bin_dir PARENT_PATH _fers_vcpkg_root_from_path)
			list(APPEND _fers_vcpkg_candidates "${_fers_vcpkg_root_from_path}")
		endif ()

		list(REMOVE_DUPLICATES _fers_vcpkg_candidates)
		foreach (_fers_vcpkg_root IN LISTS _fers_vcpkg_candidates)
			set(_fers_vcpkg_toolchain "${_fers_vcpkg_root}/scripts/buildsystems/vcpkg.cmake")
			if (EXISTS "${_fers_vcpkg_toolchain}")
				set(CMAKE_TOOLCHAIN_FILE "${_fers_vcpkg_toolchain}" CACHE FILEPATH "vcpkg toolchain")
				message(STATUS "Using vcpkg toolchain: ${CMAKE_TOOLCHAIN_FILE}")
				break()
			endif ()
		endforeach ()

		unset(_fers_vcpkg_candidates)
		unset(_fers_vcpkg_exe)
		unset(_fers_vcpkg_bin_dir)
		unset(_fers_vcpkg_root_from_path)
		unset(_fers_vcpkg_root)
		unset(_fers_vcpkg_toolchain)
	endif ()

	if (NOT DEFINED CMAKE_TOOLCHAIN_FILE)
		message(FATAL_ERROR
				"Could not locate vcpkg toolchain file. Set VCPKG_ROOT or CMAKE_TOOLCHAIN_FILE before configuring.")
	endif ()
endif ()

project(FERS
		VERSION ${FERS_PROJECT_VERSION}
		DESCRIPTION "The Flexible Extensible Radar Simulator"
		LANGUAGES CXX)

if (NOT FERS_DOCS_ONLY)
	# --- Project-wide Settings ---
	set(CMAKE_CXX_STANDARD 23)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
	set(CMAKE_CXX_EXTENSIONS OFF)
	set(CMAKE_VERBOSE_MAKEFILE ON)
	set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
	# CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS simplifies DLL exports on Windows.
	# macOS and Linux do not need this but it is harmless there.
	# Windows is not a currently supported target but this flag is kept to ease
	# a future port.
	if (WIN32)
		set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
	endif ()

	# --- Compiler Warnings and Flags ---
	include(FersCompilerWarnings)

	# --- Testing ---
	if (FERS_BUILD_TESTS)
		enable_testing()
	endif ()

	if (FERS_ENABLE_COVERAGE)
		if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
			message(STATUS "Enabling Code Coverage")
			add_compile_options(-fprofile-arcs -ftest-coverage -fprofile-update=atomic)
			add_link_options(--coverage)
		else ()
			message(WARNING "Code coverage is only supported for GCC and Clang")
		endif ()
	endif ()

	# --- macOS Specific Configuration ---
	if (APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
		message(STATUS "macOS with a modern toolchain (e.g., Homebrew LLVM) is recommended.")
		message(STATUS "If you are not using the system compiler, please configure with:")
		message(STATUS "  CC=/path/to/clang CXX=/path/to/clang++ cmake ..")
	endif ()

	# Add the subdirectory containing the C++ packages.
	add_subdirectory(packages)
endif ()

# --- Documentation ---
add_subdirectory(docs)
