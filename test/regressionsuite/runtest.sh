#!/bin/bash

base=`basename $1 .fersxml`
cd $base
../../../src/fers ../$1 
#> $base.log 2>&1
#find -maxdepth 1 -name "*.fersxml" -exec echo Diff results for: '{}' \; -exec diff -u 'correct/{}' '{}' \;