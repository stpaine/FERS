// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import {
    formatNumberFieldValue,
    resolveNumberFieldBlur,
    resolveTextFieldBlur,
} from './InspectorControls';
import {
    createDechirpReference,
    createFmcwModeConfig,
    DECHIRP_MODE_OPTIONS,
    DECHIRP_REFERENCE_SOURCE_OPTIONS,
    getAvailableDechirpReferenceSourceOptions,
    getCompatibleWaveforms,
    getFmcwEmitterNames,
    getFmcwWaveformNames,
    getPulsedRadarFieldLabels,
    RADAR_MODE_OPTIONS,
    resolveWaveformSelectValue,
    shouldClearWaveformForRadarType,
} from './PlatformComponentInspector';
import {
    ensureCubicPositionWaypoints,
    ensureCubicRotationWaypoints,
} from './PlatformInspector';
import {
    createWaveformForType,
    getVisibleWaveformFieldLabels,
    WAVEFORM_TYPE_OPTIONS,
} from './WaveformInspector';

describe('InspectorControls blur resolution', () => {
    test('reverts required numeric fields when blurred empty', () => {
        expect(resolveNumberFieldBlur('', 42, 'revert')).toEqual({
            action: 'revert',
            error: 'Value required.',
            nextDraft: '42',
        });
    });

    test('commits null for unsettable numeric fields when blurred empty', () => {
        expect(resolveNumberFieldBlur('', 42, 'null')).toEqual({
            action: 'commit',
            nextDraft: '',
            value: null,
        });
    });

    test('reverts timing noise-entry drafts instead of producing empty commits', () => {
        expect(resolveNumberFieldBlur('', 0.5, 'revert')).toEqual({
            action: 'revert',
            error: 'Value required.',
            nextDraft: '0.5',
        });
    });

    test('rejects invalid numeric drafts on blur', () => {
        expect(resolveNumberFieldBlur('-', 7, 'revert')).toEqual({
            action: 'revert',
            error: 'Invalid number.',
            nextDraft: '7',
        });
    });

    test('parses valid numeric drafts on blur', () => {
        expect(resolveNumberFieldBlur('1e3', 7, 'revert')).toEqual({
            action: 'commit',
            nextDraft: '1000',
            value: 1000,
        });
    });

    test('reverts required text fields when blurred empty', () => {
        expect(resolveTextFieldBlur('   ', 'Timing A', false)).toEqual({
            action: 'revert',
            error: 'Value required.',
            nextDraft: 'Timing A',
        });
    });

    test('commits non-empty text drafts', () => {
        expect(
            resolveTextFieldBlur('Updated Timing', 'Timing A', false)
        ).toEqual({
            action: 'commit',
            nextDraft: 'Updated Timing',
            value: 'Updated Timing',
        });
    });

    test('formats null numeric values as an empty draft', () => {
        expect(formatNumberFieldValue(null)).toBe('');
    });
});

describe('Waveform inspector authoring options', () => {
    test('offers pulse file, CW, and FMCW waveform types', () => {
        expect(WAVEFORM_TYPE_OPTIONS).toEqual([
            { value: 'pulsed_from_file', label: 'Pulse File' },
            { value: 'cw', label: 'CW' },
            { value: 'fmcw_linear_chirp', label: 'FMCW Linear Chirp' },
            { value: 'fmcw_triangle', label: 'FMCW Triangle' },
        ]);
    });

    test('creates FMCW defaults while preserving common waveform fields', () => {
        const fmcwWaveform = createWaveformForType(
            {
                id: '1',
                type: 'Waveform',
                name: 'Search chirp',
                waveformType: 'cw',
                power: 250,
                carrier_frequency: 9.6e9,
                chirp_count: 4.8,
                start_frequency_offset: 125,
            },
            'fmcw_linear_chirp'
        );

        expect(fmcwWaveform).toMatchObject({
            id: '1',
            type: 'Waveform',
            name: 'Search chirp',
            waveformType: 'fmcw_linear_chirp',
            direction: 'up',
            power: 250,
            carrier_frequency: 9.6e9,
            chirp_bandwidth: 4e3,
            chirp_duration: 1e-3,
            chirp_period: 1e-3,
            start_frequency_offset: 125,
            chirp_count: 4,
        });
    });

    test('reports FMCW chirp fields only for FMCW waveforms', () => {
        expect(getVisibleWaveformFieldLabels('fmcw_linear_chirp')).toEqual([
            'Direction',
            'Chirp Bandwidth (Hz)',
            'Chirp Duration (s)',
            'Chirp Period (s)',
            'Start Frequency Offset (Hz)',
            'Chirp Count',
        ]);
        expect(getVisibleWaveformFieldLabels('fmcw_triangle')).toEqual([
            'Chirp Bandwidth (Hz)',
            'Chirp Duration (s)',
            'Start Frequency Offset (Hz)',
            'Triangle Count',
        ]);
        expect(getVisibleWaveformFieldLabels('cw')).toEqual([]);
    });
});

describe('Platform inspector authoring helpers', () => {
    test('adds a second position waypoint before selecting cubic interpolation', () => {
        const waypoints = [
            {
                id: 'position-1',
                x: 1,
                y: 2,
                altitude: 3,
                time: 4,
            },
        ];

        const next = ensureCubicPositionWaypoints(waypoints);

        expect(next).toHaveLength(2);
        expect(next[0]).toEqual(waypoints[0]);
        expect(next[1]).toMatchObject({
            x: 1,
            y: 2,
            altitude: 3,
            time: 5,
        });
        expect(next[1].id).not.toBe('position-1');
    });

    test('adds a second rotation waypoint before selecting cubic interpolation', () => {
        const waypoints = [
            {
                id: 'rotation-1',
                azimuth: 10,
                elevation: 20,
                time: 4,
            },
        ];

        const next = ensureCubicRotationWaypoints(waypoints);

        expect(next).toHaveLength(2);
        expect(next[0]).toEqual(waypoints[0]);
        expect(next[1]).toMatchObject({
            azimuth: 10,
            elevation: 20,
            time: 5,
        });
        expect(next[1].id).not.toBe('rotation-1');
    });
});

describe('Platform component inspector waveform compatibility', () => {
    const waveforms = [
        { id: '1', name: 'Pulse', waveformType: 'pulsed_from_file' },
        { id: '2', name: 'Tone', waveformType: 'cw' },
        { id: '3', name: 'Chirp', waveformType: 'fmcw_linear_chirp' },
        { id: '4', name: 'Triangle', waveformType: 'fmcw_triangle' },
    ];

    test('offers FMCW radar mode and hides pulsed-only fields for FMCW/CW', () => {
        expect(RADAR_MODE_OPTIONS).toEqual([
            { value: 'pulsed', label: 'Pulsed' },
            { value: 'cw', label: 'CW' },
            { value: 'fmcw', label: 'FMCW' },
        ]);

        expect(getPulsedRadarFieldLabels('pulsed')).toEqual([
            'PRF (Hz)',
            'Window Skip (s)',
            'Window Length (s)',
        ]);
        expect(getPulsedRadarFieldLabels('cw')).toEqual([]);
        expect(getPulsedRadarFieldLabels('fmcw')).toEqual([]);
    });

    test('filters waveform dropdown choices by radar mode', () => {
        expect(getCompatibleWaveforms(waveforms, 'pulsed')).toEqual([
            waveforms[0],
        ]);
        expect(getCompatibleWaveforms(waveforms, 'cw')).toEqual([waveforms[1]]);
        expect(getCompatibleWaveforms(waveforms, 'fmcw')).toEqual([
            waveforms[2],
            waveforms[3],
        ]);
    });

    test('clears and blocks incompatible waveform selection values', () => {
        expect(shouldClearWaveformForRadarType('1', waveforms, 'fmcw')).toBe(
            true
        );
        expect(shouldClearWaveformForRadarType('3', waveforms, 'fmcw')).toBe(
            false
        );
        expect(shouldClearWaveformForRadarType('4', waveforms, 'fmcw')).toBe(
            false
        );
        expect(resolveWaveformSelectValue('1', waveforms, 'fmcw')).toBe('');
        expect(resolveWaveformSelectValue('3', waveforms, 'fmcw')).toBe('3');
        expect(resolveWaveformSelectValue('4', waveforms, 'fmcw')).toBe('4');
    });

    test('offers dechirp modes and source options for receiver FMCW mode', () => {
        expect(DECHIRP_MODE_OPTIONS).toEqual([
            { value: 'none', label: 'None' },
            { value: 'physical', label: 'Physical' },
            { value: 'ideal', label: 'Ideal' },
        ]);
        expect(DECHIRP_REFERENCE_SOURCE_OPTIONS).toEqual([
            { value: 'attached', label: 'Attached' },
            { value: 'transmitter', label: 'Transmitter' },
            { value: 'custom', label: 'Custom Waveform' },
        ]);
        expect(getAvailableDechirpReferenceSourceOptions('monostatic')).toEqual(
            [...DECHIRP_REFERENCE_SOURCE_OPTIONS]
        );
        expect(getAvailableDechirpReferenceSourceOptions('receiver')).toEqual([
            { value: 'transmitter', label: 'Transmitter' },
            { value: 'custom', label: 'Custom Waveform' },
        ]);
        expect(
            getAvailableDechirpReferenceSourceOptions('receiver', 'attached')
        ).toEqual([...DECHIRP_REFERENCE_SOURCE_OPTIONS]);
    });

    test('builds backend-valid dechirp config shapes', () => {
        expect(createFmcwModeConfig('none')).toEqual({});
        expect(createFmcwModeConfig('physical')).toEqual({
            dechirp_mode: 'physical',
            dechirp_reference: { source: 'attached' },
        });
        expect(
            createFmcwModeConfig('ideal', {
                dechirp_mode: 'physical',
                dechirp_reference: {
                    source: 'transmitter',
                    transmitter_name: 'TX A',
                },
            })
        ).toEqual({
            dechirp_mode: 'ideal',
            dechirp_reference: {
                source: 'transmitter',
                transmitter_name: 'TX A',
            },
        });
        expect(
            createDechirpReference('transmitter', {
                source: 'custom',
                waveform_name: 'Bad carryover',
                transmitter_name: 'TX B',
            })
        ).toEqual({
            source: 'transmitter',
            transmitter_name: 'TX B',
        });
        expect(
            createDechirpReference('custom', {
                source: 'transmitter',
                transmitter_name: 'Bad carryover',
                waveform_name: 'LO Waveform',
            })
        ).toEqual({
            source: 'custom',
            waveform_name: 'LO Waveform',
        });
    });

    test('lists only FMCW waveforms and FMCW emitters for dechirp references', () => {
        const fullWaveforms = [
            {
                id: '1',
                type: 'Waveform' as const,
                name: 'Pulse',
                waveformType: 'pulsed_from_file' as const,
                power: 1,
                carrier_frequency: 1,
                filename: 'pulse.h5',
            },
            {
                id: '2',
                type: 'Waveform' as const,
                name: 'FMCW LO',
                waveformType: 'fmcw_linear_chirp' as const,
                direction: 'up' as const,
                power: 1,
                carrier_frequency: 1,
                chirp_bandwidth: 10,
                chirp_duration: 1,
                chirp_period: 1,
                start_frequency_offset: 0,
                chirp_count: null,
            },
        ];
        const platforms = [
            {
                id: 'p1',
                type: 'Platform' as const,
                name: 'Platform',
                motionPath: {
                    interpolation: 'static' as const,
                    waypoints: [
                        {
                            id: 'wp1',
                            x: 0,
                            y: 0,
                            altitude: 0,
                            time: 0,
                        },
                    ],
                },
                rotation: {
                    type: 'fixed' as const,
                    startAzimuth: 0,
                    startElevation: 0,
                    azimuthRate: 0,
                    elevationRate: 0,
                },
                components: [
                    {
                        id: 'tx1',
                        type: 'transmitter' as const,
                        name: 'FMCW TX',
                        radarType: 'fmcw' as const,
                        prf: null,
                        antennaId: null,
                        waveformId: '2',
                        timingId: null,
                        schedule: [],
                    },
                    {
                        id: 'tx2',
                        type: 'transmitter' as const,
                        name: 'Pulse TX',
                        radarType: 'pulsed' as const,
                        prf: 1000,
                        antennaId: null,
                        waveformId: '1',
                        timingId: null,
                        schedule: [],
                    },
                ],
            },
        ];

        expect(getFmcwWaveformNames(fullWaveforms)).toEqual(['FMCW LO']);
        expect(getFmcwEmitterNames(platforms, fullWaveforms)).toEqual([
            'FMCW TX',
        ]);
    });
});
