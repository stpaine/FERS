// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { generateSimId } from './idUtils';
import { Antenna, GlobalParameters, Platform, Timing, Waveform } from './types';

type WaveformType = Waveform['waveformType'];
type WaveformDefaults<T extends WaveformType = WaveformType> =
    T extends WaveformType
        ? Omit<Extract<Waveform, { waveformType: T }>, 'id' | 'name'>
        : never;

export const defaultGlobalParameters: GlobalParameters = {
    id: 'global-parameters',
    type: 'GlobalParameters',
    rotationAngleUnit: 'deg',
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

export const defaultWaveform: WaveformDefaults<'pulsed_from_file'> = {
    type: 'Waveform',
    waveformType: 'pulsed_from_file',
    power: 1000,
    carrier_frequency: 1e9,
    filename: '',
};

export function createWaveformForType(
    waveformType: 'pulsed_from_file'
): WaveformDefaults<'pulsed_from_file'>;
export function createWaveformForType(
    waveformType: 'cw'
): WaveformDefaults<'cw'>;
export function createWaveformForType(
    waveformType: 'fmcw_linear_chirp'
): WaveformDefaults<'fmcw_linear_chirp'>;
export function createWaveformForType(
    waveformType: 'fmcw_triangle'
): WaveformDefaults<'fmcw_triangle'>;
export function createWaveformForType(
    waveformType: WaveformType
): WaveformDefaults {
    const common = {
        type: 'Waveform' as const,
        power: defaultWaveform.power,
        carrier_frequency: defaultWaveform.carrier_frequency,
    };

    switch (waveformType) {
        case 'pulsed_from_file':
            return {
                ...common,
                waveformType,
                filename: '',
            };
        case 'cw':
            return {
                ...common,
                waveformType,
            };
        case 'fmcw_linear_chirp':
            return {
                ...common,
                waveformType,
                direction: 'up',
                chirp_bandwidth: 4e3,
                chirp_duration: 1e-3,
                chirp_period: 1e-3,
                start_frequency_offset: 0,
                chirp_count: null,
            };
        case 'fmcw_triangle':
            return {
                ...common,
                waveformType,
                chirp_bandwidth: 4e3,
                chirp_duration: 1e-3,
                start_frequency_offset: 0,
                triangle_count: null,
            };
    }
}

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
