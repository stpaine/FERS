# This module defines project-wide options for FERS.

# Set a default build type for single-configuration generators if not specified.
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type" FORCE)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
	message(STATUS "No build type specified, defaulting to Release.")
endif ()

# --- Component Build Options ---
# Determine default 'ON'/'OFF' based on whether FERS is the top-level project.
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
	set(FERS_EXTRAS_DEFAULT ON)
else ()
	set(FERS_EXTRAS_DEFAULT OFF)
endif ()

option(FERS_BUILD_SHARED_LIBS "Build the shared (.so/.dll) library" ON)
option(FERS_BUILD_STATIC_LIBS "Build the static (.a) library" ON)
option(FERS_BUILD_DOCS "Enable building Doxygen documentation" OFF)
option(FERS_DOCS_ONLY "Configure only for documentation generation" OFF)
option(FERS_BUILD_TESTS "Build unit tests" ON)
option(FERS_ENABLE_COVERAGE "Enable code coverage generation" OFF)

# Ensure at least one library type is selected.
if (NOT FERS_BUILD_SHARED_LIBS AND NOT FERS_BUILD_STATIC_LIBS)
	message(FATAL_ERROR "Either FERS_BUILD_SHARED_LIBS or FERS_BUILD_STATIC_LIBS must be enabled.")
endif ()
