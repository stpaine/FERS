#!/bin/bash

find -name "test*.fersxml" -exec ./runtest.sh '{}' \;
