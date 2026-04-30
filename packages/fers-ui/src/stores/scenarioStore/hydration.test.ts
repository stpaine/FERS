// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import { defaultGlobalParameters } from './defaults';
import { buildHydratedScenarioState, parseScenarioData } from './hydration';
import { ScenarioData, ScenarioState } from './types';

function createScenarioData(): ScenarioData {
    return {
        globalParameters: {
            ...defaultGlobalParameters,
            start: 5,
            end: 20,
        },
        waveforms: [],
        timings: [],
        antennas: [],
        platforms: [
            {
                id: 'platform-1',
                type: 'Platform',
                name: 'Platform 1',
                motionPath: {
                    interpolation: 'static',
                    waypoints: [
                        {
                            id: 'waypoint-1',
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
                        id: 'component-1',
                        type: 'transmitter',
                        name: 'TX 1',
                        radarType: 'pulsed',
                        prf: 1000,
                        antennaId: null,
                        waveformId: null,
                        timingId: null,
                        schedule: [],
                    },
                ],
            },
        ],
    };
}

function createState(overrides: Partial<ScenarioState> = {}): ScenarioState {
    return {
        globalParameters: defaultGlobalParameters,
        waveforms: [],
        timings: [],
        antennas: [],
        platforms: [],
        selectedItemId: null,
        selectedComponentId: null,
        isDirty: false,
        isPlaying: false,
        currentTime: 0,
        targetPlaybackDuration: null,
        isBackendSyncing: false,
        backendVersion: 0,
        scenarioFilePath: null,
        outputDirectory: null,
        antennaPreviewErrors: {},
        notificationSnackbar: {
            open: false,
            message: '',
            severity: 'error',
        },
        notificationQueue: [],
        viewControlAction: { type: null, timestamp: 0 },
        visibility: {
            showAxes: true,
            showPatterns: true,
            showBoresights: true,
            showLinks: true,
            showLinkLabels: true,
            showLinkMonostatic: true,
            showLinkIlluminator: true,
            showLinkScattered: true,
            showLinkDirect: true,
            showVelocities: true,
            showPlatforms: true,
            showPlatformLabels: true,
            showMotionPaths: true,
        },
        ...overrides,
    };
}

describe('buildHydratedScenarioState', () => {
    test('preserves selection and clamps current time during backend recovery', () => {
        const hydrated = buildHydratedScenarioState(
            createState({
                selectedItemId: 'platform-1',
                selectedComponentId: 'component-1',
                currentTime: 99,
            }),
            createScenarioData(),
            {
                isDirty: false,
                preserveSelection: true,
                preserveCurrentTime: true,
            }
        );

        expect(hydrated.selectedItemId).toBe('platform-1');
        expect(hydrated.selectedComponentId).toBe('component-1');
        expect(hydrated.currentTime).toBe(20);
        expect(hydrated.isDirty).toBe(false);
    });

    test('clears stale selection and resets time for import-style loads', () => {
        const hydrated = buildHydratedScenarioState(
            createState({
                selectedItemId: 'missing-platform',
                selectedComponentId: 'missing-component',
                currentTime: 12,
            }),
            createScenarioData(),
            {
                isDirty: true,
            }
        );

        expect(hydrated.selectedItemId).toBeNull();
        expect(hydrated.selectedComponentId).toBeNull();
        expect(hydrated.currentTime).toBe(5);
        expect(hydrated.isDirty).toBe(true);
    });
});

describe('parseScenarioData FMCW hydration', () => {
    test('hydrates FMCW waveform and component mode from backend data', () => {
        const scenario = parseScenarioData({
            name: 'FMCW Scenario',
            parameters: {},
            waveforms: [
                {
                    id: 1,
                    name: 'FMCW Down Chirp',
                    power: 50,
                    carrier_frequency: 10e9,
                    fmcw_linear_chirp: {
                        direction: 'down',
                        chirp_bandwidth: 20e6,
                        chirp_duration: 250e-6,
                        chirp_period: 300e-6,
                        start_frequency_offset: -10e6,
                        chirp_count: 8,
                    },
                },
                {
                    id: 4,
                    name: 'FMCW Triangle',
                    power: 50,
                    carrier_frequency: 10e9,
                    fmcw_triangle: {
                        chirp_bandwidth: 20e6,
                        chirp_duration: 250e-6,
                        start_frequency_offset: -10e6,
                        triangle_count: 8,
                    },
                },
            ],
            timings: [],
            antennas: [],
            platforms: [
                {
                    id: 2,
                    name: 'Radar Platform',
                    motionpath: {
                        interpolation: 'static',
                        positionwaypoints: [
                            { x: 0, y: 0, altitude: 0, time: 0 },
                        ],
                    },
                    fixedrotation: {
                        startazimuth: 0,
                        startelevation: 0,
                        azimuthrate: 0,
                        elevationrate: 0,
                    },
                    components: [
                        {
                            receiver: {
                                id: 3,
                                name: 'FMCW RX',
                                fmcw_mode: {
                                    dechirp_mode: 'ideal',
                                    dechirp_reference: {
                                        source: 'custom',
                                        waveform_name: 'FMCW Down Chirp',
                                    },
                                },
                            },
                        },
                    ],
                },
            ],
        });

        expect(scenario).not.toBeNull();
        expect(scenario?.waveforms[0]).toMatchObject({
            id: '1',
            waveformType: 'fmcw_linear_chirp',
            direction: 'down',
            chirp_bandwidth: 20e6,
            chirp_duration: 250e-6,
            chirp_period: 300e-6,
            start_frequency_offset: -10e6,
            chirp_count: 8,
        });
        expect(scenario?.waveforms[1]).toMatchObject({
            id: '4',
            waveformType: 'fmcw_triangle',
            chirp_bandwidth: 20e6,
            chirp_duration: 250e-6,
            start_frequency_offset: -10e6,
            triangle_count: 8,
        });
        expect(scenario?.platforms[0].components[0]).toMatchObject({
            id: '3',
            type: 'receiver',
            radarType: 'fmcw',
            fmcwModeConfig: {
                dechirp_mode: 'ideal',
                dechirp_reference: {
                    source: 'custom',
                    waveform_name: 'FMCW Down Chirp',
                },
            },
        });
    });

    test('hydrates receiver and monostatic dechirp references by backend names', () => {
        const scenario = parseScenarioData({
            name: 'FMCW Dechirp Scenario',
            parameters: {},
            waveforms: [
                {
                    id: 10,
                    name: 'FMCW LO',
                    power: 50,
                    carrier_frequency: 10e9,
                    fmcw_linear_chirp: {
                        direction: 'up',
                        chirp_bandwidth: 20e6,
                        chirp_duration: 250e-6,
                        chirp_period: 300e-6,
                    },
                },
            ],
            timings: [],
            antennas: [],
            platforms: [
                {
                    id: 20,
                    name: 'Radar Platform',
                    motionpath: {
                        interpolation: 'static',
                        positionwaypoints: [
                            { x: 0, y: 0, altitude: 0, time: 0 },
                        ],
                    },
                    fixedrotation: {
                        startazimuth: 0,
                        startelevation: 0,
                        azimuthrate: 0,
                        elevationrate: 0,
                    },
                    components: [
                        {
                            transmitter: {
                                id: 21,
                                name: 'Reference TX',
                                waveform: 10,
                                fmcw_mode: {},
                            },
                        },
                        {
                            receiver: {
                                id: 22,
                                name: 'FMCW RX',
                                fmcw_mode: {
                                    dechirp_mode: 'physical',
                                    dechirp_reference: {
                                        source: 'transmitter',
                                        transmitter_name: 'Reference TX',
                                    },
                                },
                            },
                        },
                        {
                            monostatic: {
                                id: 23,
                                tx_id: 24,
                                rx_id: 25,
                                name: 'FMCW Mono',
                                waveform: 10,
                                fmcw_mode: {
                                    dechirp_mode: 'ideal',
                                    dechirp_reference: {
                                        source: 'attached',
                                    },
                                },
                            },
                        },
                    ],
                },
            ],
        });

        expect(scenario).not.toBeNull();
        expect(scenario?.platforms[0].components[1]).toMatchObject({
            id: '22',
            type: 'receiver',
            radarType: 'fmcw',
            fmcwModeConfig: {
                dechirp_mode: 'physical',
                dechirp_reference: {
                    source: 'transmitter',
                    transmitter_name: 'Reference TX',
                },
            },
        });
        expect(scenario?.platforms[0].components[2]).toMatchObject({
            id: '23',
            type: 'monostatic',
            txId: '24',
            rxId: '25',
            radarType: 'fmcw',
            fmcwModeConfig: {
                dechirp_mode: 'ideal',
                dechirp_reference: {
                    source: 'attached',
                },
            },
        });
    });
});
