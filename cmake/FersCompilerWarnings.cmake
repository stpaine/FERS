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
						   -fno-finite-math-only
	)

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
