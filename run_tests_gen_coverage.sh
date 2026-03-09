#!/bin/bash

#rm -rf build_cov/
mkdir -p build_cov

cmake -B build_cov -S . -DFERS_BUILD_TESTS=ON -DFERS_ENABLE_COVERAGE=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build_cov -j12

cd build_cov || exit
ctest --output-on-failure -j12

lcov --capture --directory . --output-file coverage.info --ignore-errors source
lcov --remove coverage.info '/usr/*' '*/_deps/*' '*/third_party/*' '*/tests/*' --output-file coverage_filtered.info
genhtml coverage_filtered.info --output-directory coverage_report --ignore-errors source
