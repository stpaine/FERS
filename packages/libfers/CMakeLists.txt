project(libfers)

include(FetchContent)

# Helper function to apply common settings to a library target
function(fers_configure_library_target target_name)
	# Apply project-wide compiler warnings and flags
	apply_fers_warnings(${target_name})

	# --- Source and Header Files ---
	target_sources(${target_name} PRIVATE
		${LIBFERS_SOURCES}
		${LIBFERS_PUBLIC_HEADERS}
		${LIBFERS_PRIVATE_HEADERS}
	)

	# --- Include Directories ---
	target_include_directories(${target_name} PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
		$<INSTALL_INTERFACE:include>
	)
	target_include_directories(${target_name} PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
	)

	# --- Schema Headers ---
	add_dependencies(${target_name} generate_schemas)
	target_include_directories(${target_name} PRIVATE
		$<BUILD_INTERFACE:${GENERATED_HEADERS_DIR}>
	)

	target_include_directories(${target_name} PRIVATE
		$<TARGET_PROPERTY:nlohmann_json::nlohmann_json,INTERFACE_INCLUDE_DIRECTORIES>
	)

	# --- Third-Party Dependencies ---
	target_link_libraries(${target_name} PRIVATE
		HighFive::HighFive
		GeographicLib::GeographicLib
		LibXml2::LibXml2
		${HDF5_C_LIBRARIES}
		${HDF5_HL_LIBRARIES}
	)
	target_compile_definitions(${target_name} PRIVATE HAVE_LIBHDF5)

	set_target_properties(${target_name} PROPERTIES
		SOVERSION ${PROJECT_VERSION_MAJOR}
		VERSION ${PROJECT_VERSION}
	)
endfunction()


# --- Gather Source and Header Files ---
include(include/libfers/CMakeLists.txt)
include(src/CMakeLists.txt)

# --- Schema Headers Generation ---
set(SCHEMA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../schema)
set(GENERATED_HEADERS_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated_headers)

add_custom_command(
	OUTPUT ${GENERATED_HEADERS_DIR}/fers_xml_dtd.h ${GENERATED_HEADERS_DIR}/fers_xml_xsd.h
	COMMAND ${CMAKE_COMMAND} -E make_directory ${GENERATED_HEADERS_DIR}
	COMMAND ${CMAKE_COMMAND} -E copy ${SCHEMA_DIR}/fers-xml.dtd ${CMAKE_CURRENT_BINARY_DIR}/fers-xml.dtd
	COMMAND ${CMAKE_COMMAND} -E copy ${SCHEMA_DIR}/fers-xml.xsd ${CMAKE_CURRENT_BINARY_DIR}/fers-xml.xsd
	COMMAND xxd -i fers-xml.dtd > ${GENERATED_HEADERS_DIR}/fers_xml_dtd.h
	COMMAND xxd -i fers-xml.xsd > ${GENERATED_HEADERS_DIR}/fers_xml_xsd.h
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	DEPENDS ${SCHEMA_DIR}/fers-xml.dtd ${SCHEMA_DIR}/fers-xml.xsd
	COMMENT "Converting XML schema files to C headers"
)
add_custom_target(generate_schemas ALL DEPENDS
	${GENERATED_HEADERS_DIR}/fers_xml_dtd.h
	${GENERATED_HEADERS_DIR}/fers_xml_xsd.h
)

# --- Third-Party Dependencies ---
add_subdirectory(third_party/highfive EXCLUDE_FROM_ALL)
add_subdirectory(third_party/json EXCLUDE_FROM_ALL)

# We configure options to disable tools, tests, and docs to speed up build
set(BUILD_DOCUMENTATION OFF CACHE BOOL "Disable GeographicLib docs" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "Disable GeographicLib tests" FORCE)
# Setting BINDIR to empty/OFF prevents GeographicLib from building CLI tools
set(BINDIR "" CACHE STRING "Disable GeographicLib tools" FORCE)
# Ensure it builds as a static library if you prefer, or let it follow BUILD_SHARED_LIBS
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build GeographicLib as static" FORCE)

FetchContent_Declare(
		GeographicLib
		GIT_REPOSITORY https://github.com/geographiclib/geographiclib.git
		GIT_TAG        r2.7
)
FetchContent_MakeAvailable(GeographicLib)

find_package(LibXml2 REQUIRED)
find_package(HDF5 COMPONENTS C HL REQUIRED)

# --- Create Library Targets ---
if (FERS_BUILD_STATIC_LIBS)
	add_library(fers_static STATIC)
	set_target_properties(fers_static PROPERTIES OUTPUT_NAME "fers")
	# Ensure the static library is built with Position-Independent Code,
	# which is required when linking it into a shared library (like the Rust cdylib).
	set_target_properties(fers_static PROPERTIES POSITION_INDEPENDENT_CODE ON)
	fers_configure_library_target(fers_static)
	add_library(fers::fers_static ALIAS fers_static)
endif ()

if (FERS_BUILD_SHARED_LIBS)
	add_library(fers_shared SHARED)
	set_target_properties(fers_shared PROPERTIES OUTPUT_NAME "fers")
	fers_configure_library_target(fers_shared)
	add_library(fers::fers_shared ALIAS fers_shared)
endif ()

# Create a generic 'libfers' INTERFACE target that downstream internal targets
# (like fers-cli) can link to. It prefers the shared library if available.
add_library(libfers INTERFACE)
if (FERS_BUILD_SHARED_LIBS)
	target_link_libraries(libfers INTERFACE fers::fers_shared)
else ()
	target_link_libraries(libfers INTERFACE fers::fers_static)
endif ()
add_library(fers::fers ALIAS libfers)

# Group sources for IDEs
source_group("Source Files" FILES ${LIBFERS_SOURCES})
source_group("Header Files\\Public" FILES ${LIBFERS_PUBLIC_HEADERS})
source_group("Header Files\\Private" FILES ${LIBFERS_PRIVATE_HEADERS})

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/third_party/json/include/nlohmann
	DESTINATION include
)

# --- Installation ---
if (FERS_BUILD_STATIC_LIBS)
	install(TARGETS fers_static
		EXPORT libfers-targets
		ARCHIVE DESTINATION lib
	)
endif ()

if (FERS_BUILD_SHARED_LIBS)
	install(TARGETS fers_shared
		EXPORT libfers-targets
		LIBRARY DESTINATION lib
		RUNTIME DESTINATION bin
		ARCHIVE DESTINATION lib
	)
endif ()

install(FILES ${LIBFERS_PUBLIC_HEADERS}
	DESTINATION include/libfers
)

install(EXPORT libfers-targets
	FILE libfers-config.cmake
	NAMESPACE fers::
	DESTINATION lib/cmake/libfers
)
