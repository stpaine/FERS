#!/bin/bash

cmake --preset=coverage
cmake --build --preset coverage
ctest --preset coverage -j12

cd build/coverage || exit 1

lcov --capture --directory . --output-file coverage.info --ignore-errors source
lcov --remove coverage.info '/usr/*' '*/_deps/*' '*/vcpkg_installed/*' '*/tests/*' --output-file coverage_filtered.info
genhtml coverage_filtered.info --output-directory coverage_report --ignore-errors source
