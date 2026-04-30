// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import { PlatformComponentSchema, WaveformSchema } from '../scenarioSchema';
import { createWaveformForType, defaultGlobalParameters } from './defaults';
import {
    serializeAntenna,
    serializeComponentInner,
    serializeGlobalParameters,
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

describe('serializeGlobalParameters', () => {
    test('fills UTM zone and hemisphere defaults before backend sync', () => {
        expect(
            serializeGlobalParameters({
                ...defaultGlobalParameters,
                coordinateSystem: { frame: 'UTM' },
            })
        ).toMatchObject({
            coordinatesystem: {
                frame: 'UTM',
                zone: 34,
                hemisphere: 'S',
            },
        });
    });
});

describe('FMCW schema', () => {
    const validWaveform: Extract<
        Waveform,
        { waveformType: 'fmcw_linear_chirp' }
    > = {
        ...createWaveformForType('fmcw_linear_chirp'),
        id: '10',
        name: 'FMCW',
    };

    test('accepts valid FMCW waveform fields', () => {
        expect(WaveformSchema.safeParse(validWaveform).success).toBe(true);
        expect(
            WaveformSchema.safeParse({
                ...createWaveformForType('fmcw_triangle'),
                id: '11',
                name: 'Triangle',
            }).success
        ).toBe(true);
    });

    test('keeps generated FMCW defaults within the backend aliasing limit', () => {
        const sweepStart = validWaveform.start_frequency_offset ?? 0;
        const sweepEnd =
            sweepStart +
            (validWaveform.direction === 'down' ? -1 : 1) *
                validWaveform.chirp_bandwidth;
        const maxBaseband = Math.max(Math.abs(sweepStart), Math.abs(sweepEnd));
        const effectiveRate =
            defaultGlobalParameters.rate *
            defaultGlobalParameters.oversample_ratio;

        expect(effectiveRate).toBeGreaterThan(maxBaseband);
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

        const invalidTriangles = [
            {
                ...createWaveformForType('fmcw_triangle'),
                id: '12',
                name: 'Bad Triangle',
                chirp_bandwidth: 0,
            },
            {
                ...createWaveformForType('fmcw_triangle'),
                id: '13',
                name: 'Bad Triangle Count',
                triangle_count: 1.5,
            },
        ];
        for (const waveform of invalidTriangles) {
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
            ...createWaveformForType('fmcw_linear_chirp'),
            id: '30',
            name: 'FMCW',
        };

        expect(serializeWaveform(waveform)).toEqual({
            id: '30',
            name: 'FMCW',
            power: 1000,
            carrier_frequency: 1e9,
            fmcw_linear_chirp: {
                direction: 'up',
                chirp_bandwidth: 4e3,
                chirp_duration: 1e-3,
                chirp_period: 1e-3,
            },
        });
    });

    test('serializes FMCW optional fields when set', () => {
        const waveform: Waveform = {
            ...createWaveformForType('fmcw_linear_chirp'),
            id: '31',
            name: 'FMCW Offset',
            start_frequency_offset: -1000,
            chirp_count: 4,
        };

        expect(serializeWaveform(waveform)).toMatchObject({
            fmcw_linear_chirp: {
                direction: 'up',
                start_frequency_offset: -1000,
                chirp_count: 4,
            },
        });
    });

    test('serializes FMCW triangle waveform payload', () => {
        const waveform: Waveform = {
            ...createWaveformForType('fmcw_triangle'),
            id: '32',
            name: 'Triangle',
            start_frequency_offset: -1000,
            triangle_count: 4,
        };

        expect(serializeWaveform(waveform)).toEqual({
            id: '32',
            name: 'Triangle',
            power: 1000,
            carrier_frequency: 1e9,
            fmcw_triangle: {
                chirp_bandwidth: 4e3,
                chirp_duration: 1e-3,
                start_frequency_offset: -1000,
                triangle_count: 4,
            },
        });
    });
});

describe('serializeComponentInner', () => {
    test('serializes empty component references as backend placeholders', () => {
        const component: PlatformComponent = {
            id: '41',
            type: 'transmitter',
            name: 'Draft Tx',
            radarType: 'pulsed',
            prf: 1000,
            antennaId: '',
            waveformId: '',
            timingId: '',
            schedule: [],
        };

        expect(serializeComponentInner(component)).toMatchObject({
            antenna: 0,
            waveform: 0,
            timing: 0,
        });
    });

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

    test('preserves receiver FMCW attached dechirp mode fields', () => {
        const component: PlatformComponent = {
            id: '42',
            type: 'receiver',
            name: 'FMCW Rx',
            radarType: 'fmcw',
            window_skip: null,
            window_length: null,
            prf: null,
            antennaId: null,
            timingId: null,
            noiseTemperature: null,
            noDirectPaths: false,
            noPropagationLoss: false,
            fmcwModeConfig: {
                dechirp_mode: 'physical',
                dechirp_reference: { source: 'attached' },
            },
            schedule: [],
        };

        expect(serializeComponentInner(component)).toMatchObject({
            fmcw_mode: {
                dechirp_mode: 'physical',
                dechirp_reference: { source: 'attached' },
            },
        });
    });

    test('preserves receiver FMCW transmitter dechirp reference name', () => {
        const component: PlatformComponent = {
            id: '43',
            type: 'receiver',
            name: 'FMCW Rx',
            radarType: 'fmcw',
            window_skip: null,
            window_length: null,
            prf: null,
            antennaId: null,
            timingId: null,
            noiseTemperature: null,
            noDirectPaths: false,
            noPropagationLoss: false,
            fmcwModeConfig: {
                dechirp_mode: 'ideal',
                dechirp_reference: {
                    source: 'transmitter',
                    transmitter_name: 'Reference TX',
                },
            },
            schedule: [],
        };

        expect(serializeComponentInner(component)).toMatchObject({
            fmcw_mode: {
                dechirp_mode: 'ideal',
                dechirp_reference: {
                    source: 'transmitter',
                    transmitter_name: 'Reference TX',
                },
            },
        });
    });

    test('preserves monostatic FMCW custom dechirp reference waveform name', () => {
        const component: PlatformComponent = {
            id: '44',
            type: 'monostatic',
            name: 'FMCW Mono',
            txId: '45',
            rxId: '46',
            radarType: 'fmcw',
            window_skip: null,
            window_length: null,
            prf: null,
            antennaId: null,
            waveformId: null,
            timingId: null,
            noiseTemperature: null,
            noDirectPaths: false,
            noPropagationLoss: false,
            fmcwModeConfig: {
                dechirp_mode: 'physical',
                dechirp_reference: {
                    source: 'custom',
                    waveform_name: 'Reference Waveform',
                },
            },
            schedule: [],
        };

        expect(serializeComponentInner(component)).toMatchObject({
            fmcw_mode: {
                dechirp_mode: 'physical',
                dechirp_reference: {
                    source: 'custom',
                    waveform_name: 'Reference Waveform',
                },
            },
        });
    });

    test('uses an isotropic backend placeholder while a file target has no filename', () => {
        const component: PlatformComponent = {
            id: '50',
            type: 'target',
            name: 'Draft Target',
            rcs_type: 'file',
            rcs_value: 3,
            rcs_model: 'constant',
        };

        expect(serializeComponentInner(component)).toEqual({
            id: '50',
            name: 'Draft Target',
            rcs: {
                type: 'isotropic',
                value: 3,
            },
        });
    });

    test('serializes file target RCS once a filename is present', () => {
        const component: PlatformComponent = {
            id: '51',
            type: 'target',
            name: 'Loaded Target',
            rcs_type: 'file',
            rcs_filename: '/tmp/target.xml',
            rcs_model: 'constant',
        };

        expect(serializeComponentInner(component)).toEqual({
            id: '51',
            name: 'Loaded Target',
            rcs: {
                type: 'file',
                filename: '/tmp/target.xml',
            },
        });
    });
});
