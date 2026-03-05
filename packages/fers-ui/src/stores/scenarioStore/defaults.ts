// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { GlobalParameters, Waveform, Timing, Antenna, Platform } from './types';
import { generateSimId } from './idUtils';

export const defaultGlobalParameters: GlobalParameters = {
    id: 'global-parameters',
    type: 'GlobalParameters',
    simulation_name: 'FERS Simulation',
    start: 0.0,
    end: 10.0,
    rate: 10000.0,
    simSamplingRate: null,
    c: 299792458.0,
    random_seed: null,
    adc_bits: 12,
    oversample_ratio: 1,
    // Default: UCT, South Africa
    origin: {
        latitude: -33.957652,
        longitude: 18.4611991,
        altitude: 111.01,
    },
    coordinateSystem: {
        frame: 'ENU',
    },
};

export const defaultWaveform: Omit<Waveform, 'id' | 'name'> = {
    type: 'Waveform',
    waveformType: 'pulsed_from_file',
    power: 1000,
    carrier_frequency: 1e9,
    filename: '',
};

export const defaultTiming: Omit<Timing, 'id' | 'name'> = {
    type: 'Timing',
    frequency: 10e6,
    freqOffset: null,
    randomFreqOffsetStdev: null,
    phaseOffset: null,
    randomPhaseOffsetStdev: null,
    noiseEntries: [],
};

export const defaultAntenna: Omit<
    Extract<Antenna, { pattern: 'isotropic' }>,
    'id' | 'name'
> = {
    type: 'Antenna',
    pattern: 'isotropic',
    efficiency: 1.0,
    meshScale: 1.0,
    design_frequency: null,
};

export const createDefaultPlatform = (): Omit<Platform, 'id' | 'name'> => ({
    type: 'Platform',
    motionPath: {
        interpolation: 'static',
        waypoints: [
            { id: generateSimId('Platform'), x: 0, y: 0, altitude: 0, time: 0 },
        ],
    },
    rotation: {
        type: 'fixed',
        startAzimuth: 0,
        startElevation: 0,
        azimuthRate: 0,
        elevationRate: 0,
    },
    components: [],
});
