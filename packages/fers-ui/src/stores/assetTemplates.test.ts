// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import {
    cloneTemplateIntoScenarioData,
    createAssetLibraryFile,
    createTemplateFromScenarioItem,
    parseAssetTemplates,
    prepareTemplatesForCatalog,
} from './assetTemplates';
import { defaultGlobalParameters } from './scenarioStore/defaults';
import type { ScenarioData } from './scenarioStore/types';

function createScenarioFixture(): ScenarioData {
    return {
        globalParameters: {
            ...defaultGlobalParameters,
            simulation_name: 'Library Fixture',
        },
        waveforms: [
            {
                id: '6001',
                type: 'Waveform',
                name: 'Pulse A',
                waveformType: 'pulsed_from_file',
                power: 100,
                carrier_frequency: 1e9,
                filename: '/tmp/pulse.h5',
            },
        ],
        timings: [
            {
                id: '7001',
                type: 'Timing',
                name: 'Timing A',
                frequency: 10e6,
                freqOffset: null,
                randomFreqOffsetStdev: null,
                phaseOffset: null,
                randomPhaseOffsetStdev: null,
                noiseEntries: [{ id: '7002', alpha: 0.2, weight: 0.5 }],
            },
        ],
        antennas: [
            {
                id: '5001',
                type: 'Antenna',
                name: 'Yagi A',
                pattern: 'sinc',
                efficiency: 0.9,
                meshScale: 1,
                design_frequency: 1e9,
                alpha: 1,
                beta: 2,
                gamma: 3,
            },
        ],
        platforms: [
            {
                id: '1001',
                type: 'Platform',
                name: 'Aircraft A',
                motionPath: {
                    interpolation: 'static',
                    waypoints: [
                        {
                            id: '1002',
                            x: 0,
                            y: 0,
                            altitude: 1000,
                            time: 0,
                        },
                    ],
                },
                rotation: {
                    type: 'path',
                    interpolation: 'static',
                    waypoints: [
                        {
                            id: '1003',
                            azimuth: 0,
                            elevation: 0,
                            time: 0,
                        },
                    ],
                },
                components: [
                    {
                        id: '2001',
                        type: 'transmitter',
                        name: 'Tx A',
                        radarType: 'pulsed',
                        prf: 1000,
                        antennaId: '5001',
                        waveformId: '6001',
                        timingId: '7001',
                        schedule: [],
                    },
                ],
                pathPoints: [
                    {
                        x: 0,
                        y: 0,
                        z: 0,
                        vx: 0,
                        vy: 0,
                        vz: 0,
                    },
                ],
                rotationPathPoints: [{ azimuth: 0, elevation: 0 }],
            },
        ],
    };
}

describe('asset templates', () => {
    test('creates top-level asset snapshots', () => {
        const scenario = createScenarioFixture();
        const template = createTemplateFromScenarioItem(scenario, '5001', {
            id: 'template-1',
            timestamp: '2026-04-13T00:00:00.000Z',
        });

        expect(template).toMatchObject({
            id: 'template-1',
            kind: 'antenna',
            name: 'Yagi A',
            payload: {
                id: '5001',
                pattern: 'sinc',
            },
        });
    });

    test('captures platform dependencies and strips runtime path caches', () => {
        const scenario = createScenarioFixture();
        const template = createTemplateFromScenarioItem(scenario, '1001', {
            id: 'template-2',
            timestamp: '2026-04-13T00:00:00.000Z',
        });

        expect(template?.kind).toBe('platform');
        if (template?.kind !== 'platform') {
            throw new Error('Expected platform template.');
        }

        expect(template.dependencies.antennas).toHaveLength(1);
        expect(template.dependencies.waveforms).toHaveLength(1);
        expect(template.dependencies.timings).toHaveLength(1);
        expect(template.payload).not.toHaveProperty('pathPoints');
        expect(template.payload).not.toHaveProperty('rotationPathPoints');
    });

    test('loads platform templates with fresh IDs and remapped dependencies', () => {
        const scenario = createScenarioFixture();
        const template = createTemplateFromScenarioItem(scenario, '1001');
        if (template?.kind !== 'platform') {
            throw new Error('Expected platform template.');
        }

        const { scenarioData, result } = cloneTemplateIntoScenarioData(
            scenario,
            template
        );
        const insertedPlatform = scenarioData.platforms.at(-1);
        const insertedComponent = insertedPlatform?.components[0];

        expect(result.warnings).toEqual([]);
        expect(insertedPlatform?.id).not.toBe('1001');
        expect(insertedPlatform?.name).toBe('Aircraft A Copy');
        expect(insertedComponent?.id).not.toBe('2001');
        expect(insertedComponent?.type).toBe('transmitter');
        if (insertedComponent?.type !== 'transmitter') {
            throw new Error('Expected transmitter component.');
        }
        expect(insertedComponent.antennaId).not.toBe('5001');
        expect(insertedComponent.waveformId).not.toBe('6001');
        expect(insertedComponent.timingId).not.toBe('7001');
        expect(
            scenarioData.antennas.some(
                (antenna) => antenna.id === insertedComponent.antennaId
            )
        ).toBe(true);
        expect(
            scenarioData.waveforms.some(
                (waveform) => waveform.id === insertedComponent.waveformId
            )
        ).toBe(true);
        expect(
            scenarioData.timings.some(
                (timing) => timing.id === insertedComponent.timingId
            )
        ).toBe(true);
    });

    test('uses template display names when loading renamed assets', () => {
        const scenario = createScenarioFixture();
        const template = createTemplateFromScenarioItem(scenario, '6001');
        if (template?.kind !== 'waveform') {
            throw new Error('Expected waveform template.');
        }

        const { scenarioData, result } = cloneTemplateIntoScenarioData(
            scenario,
            {
                ...template,
                name: 'Reusable Pulse',
            }
        );

        expect(result.insertedName).toBe('Reusable Pulse Copy');
        expect(scenarioData.waveforms.at(-1)?.name).toBe('Reusable Pulse Copy');
    });

    test('clears missing platform dependency references with warnings', () => {
        const scenario = createScenarioFixture();
        const template = createTemplateFromScenarioItem(scenario, '1001');
        if (template?.kind !== 'platform') {
            throw new Error('Expected platform template.');
        }
        const templateWithoutDependencies = {
            ...template,
            dependencies: {
                waveforms: [],
                timings: [],
                antennas: [],
            },
        };

        const { scenarioData, result } = cloneTemplateIntoScenarioData(
            scenario,
            templateWithoutDependencies
        );
        const insertedComponent = scenarioData.platforms.at(-1)?.components[0];

        expect(result.warnings).toHaveLength(3);
        expect(insertedComponent?.type).toBe('transmitter');
        if (insertedComponent?.type !== 'transmitter') {
            throw new Error('Expected transmitter component.');
        }
        expect(insertedComponent.antennaId).toBeNull();
        expect(insertedComponent.waveformId).toBeNull();
        expect(insertedComponent.timingId).toBeNull();
    });

    test('parses catalog files and assigns imported templates new catalog IDs', () => {
        const scenario = createScenarioFixture();
        const template = createTemplateFromScenarioItem(scenario, '6001', {
            id: 'template-original',
        });
        if (!template) {
            throw new Error('Expected waveform template.');
        }

        const parsed = parseAssetTemplates(createAssetLibraryFile([template]));
        const prepared = prepareTemplatesForCatalog(parsed);

        expect(parsed).toHaveLength(1);
        expect(prepared).toHaveLength(1);
        expect(prepared[0].id).not.toBe('template-original');
        expect(prepared[0].payload.id).toBe('6001');
    });
});
