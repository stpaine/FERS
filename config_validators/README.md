# Validators

Thanks Michael Altshuler, your effort did not unnoticed

## Purpose

Allows one to check a fersxml file is configured correctly and allows one to generate a kml file used for plotting on a map

## Requirements

Ensure one has run `sudo apt install libxerces-c-dev` to ensure it is availble to be incldued and linked.

## Building

Unlike FERS above, this uses cmake to simplify building

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. run `make`

## Usage

1. run `validator <file.fersxml>'
2. follow prompts