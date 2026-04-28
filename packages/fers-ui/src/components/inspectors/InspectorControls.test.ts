// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import {
    formatNumberFieldValue,
    resolveNumberFieldBlur,
    resolveTextFieldBlur,
} from './InspectorControls';
import {
    getCompatibleWaveforms,
    getPulsedRadarFieldLabels,
    RADAR_MODE_OPTIONS,
    resolveWaveformSelectValue,
    shouldClearWaveformForRadarType,
} from './PlatformComponentInspector';
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
    test('offers pulse file, CW, and FMCW linear chirp waveform types', () => {
        expect(WAVEFORM_TYPE_OPTIONS).toEqual([
            { value: 'pulsed_from_file', label: 'Pulse File' },
            { value: 'cw', label: 'CW' },
            { value: 'fmcw_linear_chirp', label: 'FMCW Linear Chirp' },
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
        expect(getVisibleWaveformFieldLabels('cw')).toEqual([]);
    });
});

describe('Platform component inspector waveform compatibility', () => {
    const waveforms = [
        { id: '1', name: 'Pulse', waveformType: 'pulsed_from_file' },
        { id: '2', name: 'Tone', waveformType: 'cw' },
        { id: '3', name: 'Chirp', waveformType: 'fmcw_linear_chirp' },
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
        ]);
    });

    test('clears and blocks incompatible waveform selection values', () => {
        expect(shouldClearWaveformForRadarType('1', waveforms, 'fmcw')).toBe(
            true
        );
        expect(shouldClearWaveformForRadarType('3', waveforms, 'fmcw')).toBe(
            false
        );
        expect(resolveWaveformSelectValue('1', waveforms, 'fmcw')).toBe('');
        expect(resolveWaveformSelectValue('3', waveforms, 'fmcw')).toBe('3');
    });
});
