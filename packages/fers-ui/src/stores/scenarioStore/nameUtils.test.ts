// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { beforeEach, describe, expect, test } from 'bun:test';
import { defaultGlobalParameters } from './defaults';
import { useScenarioStore } from './index';
import {
    resetSyncQueueForTests,
    setSyncQueueInvokerForTests,
    waitForSyncIdle,
} from './syncQueue';

type InvokeFn = typeof import('@tauri-apps/api/core').invoke;

describe('scenario object names', () => {
    beforeEach(() => {
        resetSyncQueueForTests();
        setSyncQueueInvokerForTests((async () => []) as InvokeFn);
        useScenarioStore.setState({
            globalParameters: defaultGlobalParameters,
            waveforms: [],
            timings: [],
            antennas: [],
            platforms: [],
            selectedItemId: null,
            selectedComponentId: null,
            isDirty: false,
            currentTime: defaultGlobalParameters.start,
        });
    });

    test('component creation auto-renames duplicate default names', async () => {
        const store = useScenarioStore.getState();
        store.addPlatform();
        const platform = useScenarioStore.getState().platforms[0];

        useScenarioStore
            .getState()
            .addPlatformComponent(platform.id, 'monostatic');
        useScenarioStore
            .getState()
            .addPlatformComponent(platform.id, 'monostatic');

        expect(
            useScenarioStore
                .getState()
                .platforms[0].components.map((component) => component.name)
        ).toEqual(['Platform 1 Monostatic', 'Platform 1 Monostatic Copy']);

        await waitForSyncIdle();
    });

    test('manual duplicate component rename gets a copy suffix', async () => {
        useScenarioStore.setState({
            waveforms: [
                {
                    id: '1688849860263937',
                    type: 'Waveform',
                    name: 'Shared Name',
                    waveformType: 'cw',
                    power: 1,
                    carrier_frequency: 1e9,
                },
            ],
            platforms: [
                {
                    id: '281474976710657',
                    type: 'Platform',
                    name: 'Platform 1',
                    motionPath: {
                        interpolation: 'static',
                        waypoints: [
                            {
                                id: '281474976710658',
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
                            id: '562949953421313',
                            type: 'transmitter',
                            name: 'Original Tx',
                            radarType: 'cw',
                            prf: null,
                            antennaId: null,
                            waveformId: null,
                            timingId: null,
                            schedule: [],
                        },
                    ],
                },
            ],
        });

        useScenarioStore
            .getState()
            .updateItem('281474976710657', 'components.0.name', 'Shared Name');

        expect(
            useScenarioStore.getState().platforms[0].components[0].name
        ).toBe('Shared Name Copy');

        await waitForSyncIdle();
    });

    test('default asset names avoid cross-kind conflicts', async () => {
        useScenarioStore.setState({
            platforms: [
                {
                    id: '281474976710657',
                    type: 'Platform',
                    name: 'Platform 1',
                    motionPath: {
                        interpolation: 'static',
                        waypoints: [
                            {
                                id: '281474976710658',
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
                            id: '1125899906842625',
                            type: 'target',
                            name: 'Waveform 1',
                            rcs_type: 'isotropic',
                            rcs_value: 1,
                            rcs_model: 'constant',
                        },
                    ],
                },
            ],
        });

        useScenarioStore.getState().addWaveform();

        expect(useScenarioStore.getState().waveforms[0].name).toBe(
            'Waveform 1 Copy'
        );

        await waitForSyncIdle();
    });
});
