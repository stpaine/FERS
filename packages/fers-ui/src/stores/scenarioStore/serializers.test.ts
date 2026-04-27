// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import { PlatformComponentSchema, WaveformSchema } from '../scenarioSchema';
import { createWaveformForType } from './defaults';
import {
    serializeAntenna,
    serializeComponentInner,
    serializeWaveform,
} from './serializers';
import type { Antenna, PlatformComponent, Waveform } from './types';

describe('serializeAntenna', () => {
    test('uses an isotropic backend placeholder while an H5 antenna has no filename', () => {
        const antenna: Antenna = {
            id: '1',
            type: 'Antenna',
            name: 'Draft H5',
            pattern: 'file',
            filename: '   ',
            efficiency: 0.8,
            meshScale: 1,
            design_frequency: null,
        };

        expect(serializeAntenna(antenna)).toEqual({
            id: '1',
            name: 'Draft H5',
            pattern: 'isotropic',
            efficiency: 0.8,
        });
    });

    test('preserves H5 antenna payload once a filename is present', () => {
        const antenna: Antenna = {
            id: '2',
            type: 'Antenna',
            name: 'Loaded H5',
            pattern: 'file',
            filename: '/tmp/pattern.h5',
            efficiency: 1,
            meshScale: 1,
            design_frequency: null,
        };

        expect(serializeAntenna(antenna)).toEqual({
            id: '2',
            name: 'Loaded H5',
            pattern: 'file',
            filename: '/tmp/pattern.h5',
            efficiency: 1,
        });
    });
});

describe('FMCW schema', () => {
    const validWaveform: Waveform = {
        ...createWaveformForType('fmcw_up_chirp'),
        id: '10',
        name: 'FMCW',
    };

    test('accepts valid FMCW waveform fields', () => {
        expect(WaveformSchema.safeParse(validWaveform).success).toBe(true);
    });

    test('rejects invalid FMCW waveform fields', () => {
        const invalidWaveforms = [
            { ...validWaveform, chirp_bandwidth: 0 },
            { ...validWaveform, chirp_duration: 0 },
            {
                ...validWaveform,
                chirp_duration: 2e-6,
                chirp_period: 1e-6,
            },
            { ...validWaveform, chirp_count: 0 },
            { ...validWaveform, chirp_count: 1.5 },
            {
                ...validWaveform,
                start_frequency_offset: Number.POSITIVE_INFINITY,
            },
        ];

        for (const waveform of invalidWaveforms) {
            expect(WaveformSchema.safeParse(waveform).success).toBe(false);
        }
    });

    test('accepts FMCW radar mode on radar components', () => {
        const component: PlatformComponent = {
            id: '20',
            type: 'transmitter',
            name: 'FMCW TX',
            radarType: 'fmcw',
            prf: null,
            antennaId: null,
            waveformId: null,
            timingId: null,
            schedule: [],
        };

        expect(PlatformComponentSchema.safeParse(component).success).toBe(true);
    });
});

describe('serializeWaveform', () => {
    test('serializes FMCW waveform payload and omits unset optional fields', () => {
        const waveform: Waveform = {
            ...createWaveformForType('fmcw_up_chirp'),
            id: '30',
            name: 'FMCW',
        };

        expect(serializeWaveform(waveform)).toEqual({
            id: '30',
            name: 'FMCW',
            power: 1000,
            carrier_frequency: 1e9,
            fmcw_up_chirp: {
                chirp_bandwidth: 20e6,
                chirp_duration: 250e-6,
                chirp_period: 250e-6,
            },
        });
    });

    test('serializes FMCW optional fields when set', () => {
        const waveform: Waveform = {
            ...createWaveformForType('fmcw_up_chirp'),
            id: '31',
            name: 'FMCW Offset',
            start_frequency_offset: -1000,
            chirp_count: 4,
        };

        expect(serializeWaveform(waveform)).toMatchObject({
            fmcw_up_chirp: {
                start_frequency_offset: -1000,
                chirp_count: 4,
            },
        });
    });
});

describe('serializeComponentInner', () => {
    test('serializes FMCW mode without pulsed fields', () => {
        const component: PlatformComponent = {
            id: '40',
            type: 'monostatic',
            name: 'FMCW Mono',
            txId: '41',
            rxId: '42',
            radarType: 'fmcw',
            window_skip: 1,
            window_length: 2,
            prf: 3,
            antennaId: null,
            waveformId: null,
            timingId: null,
            noiseTemperature: null,
            noDirectPaths: false,
            noPropagationLoss: false,
            schedule: [],
        };

        expect(serializeComponentInner(component)).toEqual({
            tx_id: '41',
            rx_id: '42',
            name: 'FMCW Mono',
            fmcw_mode: {},
            antenna: 0,
            waveform: 0,
            timing: 0,
            noise_temp: null,
            nodirect: false,
            nopropagationloss: false,
            schedule: [],
        });
    });
});
