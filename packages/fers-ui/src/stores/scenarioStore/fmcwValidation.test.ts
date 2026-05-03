// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import { createWaveformForType, defaultGlobalParameters } from './defaults';
import { validateFmcwScenario, validateFmcwWaveform } from './fmcwValidation';
import type { ScenarioData, Waveform } from './types';

const makeLinearWaveform = (
    overrides: Partial<Extract<Waveform, { waveformType: 'fmcw_linear_chirp' }>>
): Extract<Waveform, { waveformType: 'fmcw_linear_chirp' }> => ({
    ...createWaveformForType('fmcw_linear_chirp'),
    id: 'wave-1',
    name: 'Linear',
    ...overrides,
});

const makeTriangleWaveform = (
    overrides: Partial<Extract<Waveform, { waveformType: 'fmcw_triangle' }>>
): Extract<Waveform, { waveformType: 'fmcw_triangle' }> => ({
    ...createWaveformForType('fmcw_triangle'),
    id: 'wave-1',
    name: 'Triangle',
    ...overrides,
});

const makeScenario = (waveform: Waveform): ScenarioData => ({
    globalParameters: defaultGlobalParameters,
    waveforms: [waveform],
    timings: [],
    antennas: [],
    platforms: [
        {
            id: 'platform-1',
            type: 'Platform',
            name: 'Platform',
            motionPath: {
                interpolation: 'static',
                waypoints: [
                    {
                        id: 'wp-1',
                        x: 0,
                        y: 0,
                        altitude: 0,
                        time: 0,
                    },
                ],
            },
            rotation: {
                type: 'fixed',
                startAzimuth: 0,
                startElevation: 0,
                azimuthRate: 0,
                elevationRate: 0,
            },
            components: [
                {
                    id: 'tx-1',
                    type: 'transmitter',
                    name: 'FMCW Tx',
                    radarType: 'fmcw',
                    prf: null,
                    antennaId: null,
                    waveformId: waveform.id,
                    timingId: null,
                    schedule: [],
                },
            ],
        },
    ],
});

describe('FMCW validation', () => {
    test('matches backend baseband and RF sweep checks', () => {
        const aliasIssues = validateFmcwWaveform(
            makeLinearWaveform({ chirp_bandwidth: 20e3 }),
            defaultGlobalParameters
        );
        expect(aliasIssues.some((issue) => issue.severity === 'error')).toBe(
            true
        );

        const rfIssues = validateFmcwWaveform(
            makeLinearWaveform({
                carrier_frequency: 100,
                start_frequency_offset: -200,
            }),
            defaultGlobalParameters
        );
        expect(
            rfIssues.some((issue) => issue.message.includes('lower sweep edge'))
        ).toBe(true);
    });

    test('reports linear schedule errors and warnings', () => {
        const scenario = makeScenario(
            makeLinearWaveform({
                chirp_duration: 1e-3,
                chirp_period: 2e-3,
            })
        );
        const component = scenario.platforms[0].components[0];
        if (component.type !== 'transmitter') {
            throw new Error('Expected transmitter fixture');
        }
        component.schedule = [
            { start: 0, end: 0.5e-3 },
            { start: 1, end: 1.0015 },
        ];

        const issues = validateFmcwScenario(scenario);
        expect(issues.some((issue) => issue.severity === 'error')).toBe(true);
        expect(issues.some((issue) => issue.severity === 'warning')).toBe(true);
    });

    test('reports triangle period errors and silent-tail warnings', () => {
        const scenario = makeScenario(
            makeTriangleWaveform({ chirp_duration: 1e-3 })
        );
        const component = scenario.platforms[0].components[0];
        if (component.type !== 'transmitter') {
            throw new Error('Expected transmitter fixture');
        }
        component.schedule = [
            { start: 0, end: 1e-3 },
            { start: 1, end: 1.0025 },
        ];

        const issues = validateFmcwScenario(scenario);
        expect(
            issues.some((issue) => issue.message.includes('triangle period'))
        ).toBe(true);
        expect(issues.some((issue) => issue.message.includes('silent'))).toBe(
            true
        );
    });

    test('reports receiver dechirp mode/reference mismatches', () => {
        const scenario = makeScenario(makeLinearWaveform({}));
        scenario.platforms[0].components.push({
            id: 'rx-1',
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
                dechirp_mode: 'none',
                dechirp_reference: { source: 'attached' },
            },
            schedule: [],
        });

        expect(
            validateFmcwScenario(scenario).some((issue) =>
                issue.message.includes('while dechirp mode is none')
            )
        ).toBe(true);

        const receiver = scenario.platforms[0].components[1];
        if (receiver.type !== 'receiver') {
            throw new Error('Expected receiver fixture');
        }
        receiver.fmcwModeConfig = { dechirp_mode: 'physical' };
        expect(
            validateFmcwScenario(scenario).some((issue) =>
                issue.message.includes('does not declare a dechirp reference')
            )
        ).toBe(true);

        receiver.fmcwModeConfig = {
            dechirp_mode: 'ideal',
            dechirp_reference: { source: 'attached' },
        };
        expect(
            validateFmcwScenario(scenario).some((issue) =>
                issue.message.includes('only monostatic receivers')
            )
        ).toBe(true);
    });

    test('validates receiver IF-chain settings', () => {
        const scenario = makeScenario(makeLinearWaveform({}));
        scenario.platforms[0].components.push({
            id: 'rx-1',
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
                dechirp_reference: {
                    source: 'transmitter',
                    transmitter_name: 'FMCW Tx',
                },
                if_sample_rate: 1e6,
                if_filter_bandwidth: 4e5,
                if_filter_transition_width: 1e5,
            },
            schedule: [],
        });

        expect(
            validateFmcwScenario(scenario).some(
                (issue) => issue.field === 'fmcwModeConfig'
            )
        ).toBe(false);

        const receiver = scenario.platforms[0].components[1];
        if (receiver.type !== 'receiver') {
            throw new Error('Expected receiver fixture');
        }

        receiver.fmcwModeConfig = {
            dechirp_mode: 'none',
            if_sample_rate: 1e6,
        };
        expect(
            validateFmcwScenario(scenario).some((issue) =>
                issue.message.includes('IF-chain settings')
            )
        ).toBe(true);

        receiver.fmcwModeConfig = {
            dechirp_mode: 'physical',
            dechirp_reference: {
                source: 'transmitter',
                transmitter_name: 'FMCW Tx',
            },
            if_sample_rate: -1,
        };
        expect(
            validateFmcwScenario(scenario).some((issue) =>
                issue.message.includes('IF sample rate')
            )
        ).toBe(true);

        receiver.fmcwModeConfig = {
            dechirp_mode: 'physical',
            dechirp_reference: {
                source: 'transmitter',
                transmitter_name: 'FMCW Tx',
            },
            if_sample_rate: 1e6,
            if_filter_bandwidth: 5e5,
        };
        expect(
            validateFmcwScenario(scenario).some((issue) =>
                issue.message.includes('less than half')
            )
        ).toBe(true);
    });

    test('validates named transmitter and custom waveform dechirp references', () => {
        const waveform = makeLinearWaveform({ id: 'wave-1', name: 'FMCW LO' });
        const scenario = makeScenario(waveform);
        scenario.waveforms.push({
            id: 'wave-2',
            type: 'Waveform',
            name: 'CW Tone',
            waveformType: 'cw',
            power: 1,
            carrier_frequency: 1e9,
        });
        scenario.platforms[0].components.push({
            id: 'rx-1',
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
                dechirp_reference: {
                    source: 'transmitter',
                    transmitter_name: 'FMCW Tx',
                },
            },
            schedule: [],
        });

        expect(
            validateFmcwScenario(scenario).some(
                (issue) => issue.field === 'fmcwModeConfig'
            )
        ).toBe(false);

        const receiver = scenario.platforms[0].components[1];
        if (receiver.type !== 'receiver') {
            throw new Error('Expected receiver fixture');
        }

        receiver.fmcwModeConfig = {
            dechirp_mode: 'physical',
            dechirp_reference: {
                source: 'transmitter',
                transmitter_name: 'Missing TX',
            },
        };
        expect(
            validateFmcwScenario(scenario).some((issue) =>
                issue.message.includes('must be an FMCW transmitter')
            )
        ).toBe(true);

        receiver.fmcwModeConfig = {
            dechirp_mode: 'ideal',
            dechirp_reference: {
                source: 'custom',
                waveform_name: 'FMCW LO',
            },
        };
        expect(
            validateFmcwScenario(scenario).some(
                (issue) => issue.field === 'fmcwModeConfig'
            )
        ).toBe(false);

        receiver.fmcwModeConfig = {
            dechirp_mode: 'ideal',
            dechirp_reference: {
                source: 'custom',
                waveform_name: 'CW Tone',
            },
        };
        expect(
            validateFmcwScenario(scenario).some((issue) =>
                issue.message.includes('must be a top-level FMCW waveform')
            )
        ).toBe(true);
    });
});
