if (NOT FERS_BUILD_DOCS)
	return()
endif()

find_package(Doxygen REQUIRED COMPONENTS dot)
if (NOT Doxygen_FOUND)
	message(WARNING "Doxygen not found, documentation will not be built. Please install Doxygen and Graphviz.")
	return()
endif()

message(STATUS "Doxygen found, documentation target will be available.")

# Define paths relative to the project root
set(DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/docs")
set(DOXYGEN_AWESOME_DIR "${CMAKE_CURRENT_SOURCE_DIR}/doxygen-awesome-css")

# Configure custom header for awesome-css extensions
set(DOXYGEN_HEADER_TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/header.html.in")
set(DOXYGEN_HEADER_FILE "${DOXYGEN_OUTPUT_DIR}/header.html")
configure_file(${DOXYGEN_HEADER_TEMPLATE} ${DOXYGEN_HEADER_FILE} @ONLY)

# Define paths to theme files for Doxyfile
set(DOXYGEN_AWESOME_CSS "${DOXYGEN_AWESOME_DIR}/doxygen-awesome.css")
set(DOXYGEN_CUSTOM_CSS "${CMAKE_CURRENT_SOURCE_DIR}/custom.css")
set(DOXYGEN_AWESOME_DARKMODE_TOGGLE_JS "${DOXYGEN_AWESOME_DIR}/doxygen-awesome-darkmode-toggle.js")
set(DOXYGEN_AWESOME_FRAGMENT_COPY_JS "${DOXYGEN_AWESOME_DIR}/doxygen-awesome-fragment-copy-button.js")
set(DOXYGEN_AWESOME_PARA_LINK_JS "${DOXYGEN_AWESOME_DIR}/doxygen-awesome-paragraph-link.js")
set(DOXYGEN_AWESOME_TABS_JS "${DOXYGEN_AWESOME_DIR}/doxygen-awesome-tabs.js")
set(DOXYGEN_AWESOME_TOC_JS "${DOXYGEN_AWESOME_DIR}/doxygen-awesome-interactive-toc.js")

# Doxygen copies these files to the output directory
set(DOXYGEN_EXTRA_FILES
	${DOXYGEN_AWESOME_DARKMODE_TOGGLE_JS}
	${DOXYGEN_AWESOME_FRAGMENT_COPY_JS}
	${DOXYGEN_AWESOME_PARA_LINK_JS}
	${DOXYGEN_AWESOME_TABS_JS}
	${DOXYGEN_AWESOME_TOC_JS}
)

# Convert the CMake list to a space-separated string for Doxygen
string(REPLACE ";" " " DOXYGEN_EXTRA_FILES_STR "${DOXYGEN_EXTRA_FILES}")

# Configure Doxyfile
set(FERS_DOXYGEN_HAVE_DOT ${Doxygen_DOT_FOUND})
configure_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile"
	"${DOXYGEN_OUTPUT_DIR}/Doxyfile"
	@ONLY
)

add_custom_target(doc
	COMMAND Doxygen::doxygen Doxyfile
	WORKING_DIRECTORY "${DOXYGEN_OUTPUT_DIR}"
	COMMENT "Generating Doxygen documentation for FERS..."
	VERBATIM
)
