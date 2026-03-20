# This module defines a function to apply consistent compiler warnings and flags
# across all project targets.

function(apply_fers_warnings TARGET_NAME)
	# Find the Threads package to enable the Threads::Threads target.
	find_package(Threads REQUIRED)

	# --- Common Warnings for GCC and Clang ---
	set(FERS_GCC_CLANG_WARNINGS
		-Wall
		-Wextra
		-Wshadow
		-Wnon-virtual-dtor
		-Wold-style-cast
		-Wcast-align
		-Wunused
		-Woverloaded-virtual
		-Wpedantic
		-Wconversion
		-Wsign-conversion
		-Wnull-dereference
		-Wdouble-promotion
		-Wformat=2
		-Wmisleading-indentation
	)

	target_compile_options(${TARGET_NAME} PRIVATE ${FERS_GCC_CLANG_WARNINGS})

	# --- GCC-Specific Warnings ---
	if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		target_compile_options(${TARGET_NAME} PRIVATE
							   -Wduplicated-cond
							   -Wduplicated-branches
							   -Wlogical-op
		)
	endif ()

	# --- Common Flags ---
	target_compile_options(${TARGET_NAME} PRIVATE
						   -pthread
						   -ffast-math
						   -fno-finite-math-only
	)

	if (UNIX)
		target_compile_definitions(${TARGET_NAME} PRIVATE _REENTRANT)
	endif ()

	# --- Build-Type Specific Flags ---
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		target_compile_options(${TARGET_NAME} PRIVATE -g -O0)
	elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
		target_compile_options(${TARGET_NAME} PRIVATE -O2)
		target_link_options(${TARGET_NAME} PRIVATE -s)
	elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
		target_compile_options(${TARGET_NAME} PRIVATE -O2 -g)
	endif ()

	# --- Sanitizer Detection ---
	if (CMAKE_CXX_FLAGS MATCHES ".*-fsanitize=address.*")
		target_compile_options(${TARGET_NAME} PRIVATE -fsanitize=address -fno-omit-frame-pointer)
		target_link_options(${TARGET_NAME} PRIVATE -fsanitize=address)
	elseif (CMAKE_CXX_FLAGS MATCHES ".*-fsanitize=thread.*")
		target_compile_options(${TARGET_NAME} PRIVATE -fsanitize=thread)
		target_link_options(${TARGET_NAME} PRIVATE -fsanitize=thread)
	endif ()

	# --- Link Threads ---
	target_link_libraries(${TARGET_NAME} PRIVATE Threads::Threads)
endfunction()
