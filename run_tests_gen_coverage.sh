#!/bin/bash

cmake --preset=coverage
cmake --build --preset=coverage -j12
ctest --preset=coverage -j12

cd build/coverage || exit 1

lcov --capture --directory . --output-file coverage.info --ignore-errors source,gcov,negative
lcov --remove coverage.info '/usr/*' '*/tests/*' '*/_deps/*' '*/vcpkg_installed/*' --output-file coverage_filtered.info
genhtml coverage_filtered.info --output-directory coverage_report --ignore-errors source
