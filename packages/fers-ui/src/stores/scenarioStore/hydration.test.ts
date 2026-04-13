// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import { defaultGlobalParameters } from './defaults';
import { buildHydratedScenarioState } from './hydration';
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
