# This module defines a function to apply consistent compiler warnings and flags
# across all project targets.

function(apply_fers_warnings TARGET_NAME)
	# Find the Threads package to enable the Threads::Threads target.
	find_package(Threads REQUIRED)

	if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
		# --- Common Warnings for GCC and Clang ---
		target_compile_options(${TARGET_NAME} PRIVATE
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
							   -fno-finite-math-only
		)

		# --- GCC-Specific Warnings ---
		if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
			target_compile_options(${TARGET_NAME} PRIVATE
								   -Wduplicated-cond
								   -Wduplicated-branches
								   -Wlogical-op
			)
		endif ()
	elseif (MSVC)
		target_compile_options(${TARGET_NAME} PRIVATE
							   /W4
							   /permissive-
							   /Zc:__cplusplus
							   /Zc:preprocessor
							   /utf-8
							   /EHsc
							   /bigobj
		)
	endif ()

	# --- Sanitizer Detection ---
	if (MSVC AND CMAKE_CXX_FLAGS MATCHES ".*\\/fsanitize=address.*")
		target_compile_options(${TARGET_NAME} PRIVATE /fsanitize=address)
		target_link_options(${TARGET_NAME} PRIVATE /fsanitize=address)
	elseif (CMAKE_CXX_FLAGS MATCHES ".*-fsanitize=address.*")
		target_compile_options(${TARGET_NAME} PRIVATE -fsanitize=address -fno-omit-frame-pointer)
		target_link_options(${TARGET_NAME} PRIVATE -fsanitize=address)
	elseif (NOT MSVC AND CMAKE_CXX_FLAGS MATCHES ".*-fsanitize=thread.*")
		target_compile_options(${TARGET_NAME} PRIVATE -fsanitize=thread)
		target_link_options(${TARGET_NAME} PRIVATE -fsanitize=thread)
	endif ()

	# --- Link Threads ---
	target_link_libraries(${TARGET_NAME} PRIVATE Threads::Threads)
endfunction()
