// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { z } from 'zod';
import {
    AntennaSchema,
    PlatformSchema,
    TimingSchema,
    WaveformSchema,
} from './scenarioSchema';
import { generateSimId } from './scenarioStore/idUtils';
import {
    collectScenarioNames,
    createUniqueName,
    createUniqueScenarioCopyName,
} from './scenarioStore/nameUtils';
import type {
    Antenna,
    Platform,
    PlatformComponent,
    ScenarioData,
    Timing,
    Waveform,
} from './scenarioStore/types';

export const ASSET_LIBRARY_SCHEMA_VERSION = 1 as const;

const TemplateBaseSchema = z.object({
    schemaVersion: z.literal(ASSET_LIBRARY_SCHEMA_VERSION),
    id: z.string().min(1),
    name: z.string().min(1),
    createdAt: z.string().min(1),
    updatedAt: z.string().min(1),
    source: z
        .object({
            scenarioName: z.string().optional(),
            itemId: z.string().optional(),
        })
        .optional(),
});

const PlatformDependenciesSchema = z.object({
    waveforms: z.array(WaveformSchema),
    timings: z.array(TimingSchema),
    antennas: z.array(AntennaSchema),
});

const WaveformTemplateSchema = TemplateBaseSchema.extend({
    kind: z.literal('waveform'),
    payload: WaveformSchema,
});

const TimingTemplateSchema = TemplateBaseSchema.extend({
    kind: z.literal('timing'),
    payload: TimingSchema,
});

const AntennaTemplateSchema = TemplateBaseSchema.extend({
    kind: z.literal('antenna'),
    payload: AntennaSchema,
});

const PlatformTemplateSchema = TemplateBaseSchema.extend({
    kind: z.literal('platform'),
    payload: PlatformSchema,
    dependencies: PlatformDependenciesSchema,
});

export const AssetLibraryTemplateSchema = z.discriminatedUnion('kind', [
    WaveformTemplateSchema,
    TimingTemplateSchema,
    AntennaTemplateSchema,
    PlatformTemplateSchema,
]);

export const AssetLibraryFileSchema = z.object({
    schemaVersion: z.literal(ASSET_LIBRARY_SCHEMA_VERSION),
    templates: z.array(AssetLibraryTemplateSchema),
});

export type AssetTemplateKind = z.infer<
    typeof AssetLibraryTemplateSchema
>['kind'];
export type AssetLibraryTemplate = z.infer<typeof AssetLibraryTemplateSchema>;
export type AssetLibraryFile = z.infer<typeof AssetLibraryFileSchema>;

export interface AssetTemplateInsertionResult {
    insertedItemId: string | null;
    insertedName: string | null;
    warnings: string[];
}

type CreateTemplateOptions = {
    id?: string;
    timestamp?: string;
};

let fallbackTemplateCounter = 1;

export function createAssetTemplateId(): string {
    if (globalThis.crypto?.randomUUID) {
        return globalThis.crypto.randomUUID();
    }

    const id = `template-${Date.now()}-${fallbackTemplateCounter}`;
    fallbackTemplateCounter += 1;
    return id;
}

function deepClone<T>(value: T): T {
    return JSON.parse(JSON.stringify(value)) as T;
}

function createBaseTemplate(
    item: Waveform | Timing | Antenna | Platform,
    scenarioData: ScenarioData,
    options: CreateTemplateOptions
) {
    const timestamp = options.timestamp ?? new Date().toISOString();
    return {
        schemaVersion: ASSET_LIBRARY_SCHEMA_VERSION,
        id: options.id ?? createAssetTemplateId(),
        name: item.name,
        createdAt: timestamp,
        updatedAt: timestamp,
        source: {
            scenarioName: scenarioData.globalParameters.simulation_name,
            itemId: item.id,
        },
    };
}

function sanitizePlatform(platform: Platform): Platform {
    const cloned = deepClone(platform);
    delete cloned.pathPoints;
    delete cloned.rotationPathPoints;
    return cloned;
}

function uniqueById<T extends { id: string }>(items: T[]): T[] {
    const seen = new Set<string>();
    return items.filter((item) => {
        if (seen.has(item.id)) {
            return false;
        }
        seen.add(item.id);
        return true;
    });
}

function collectPlatformDependencies(
    platform: Platform,
    scenarioData: ScenarioData
) {
    const antennaIds = new Set<string>();
    const waveformIds = new Set<string>();
    const timingIds = new Set<string>();

    for (const component of platform.components) {
        if ('antennaId' in component && component.antennaId) {
            antennaIds.add(component.antennaId);
        }
        if ('waveformId' in component && component.waveformId) {
            waveformIds.add(component.waveformId);
        }
        if ('timingId' in component && component.timingId) {
            timingIds.add(component.timingId);
        }
    }

    return {
        waveforms: uniqueById(
            scenarioData.waveforms.filter((waveform) =>
                waveformIds.has(waveform.id)
            )
        ).map(deepClone),
        timings: uniqueById(
            scenarioData.timings.filter((timing) => timingIds.has(timing.id))
        ).map(deepClone),
        antennas: uniqueById(
            scenarioData.antennas.filter((antenna) =>
                antennaIds.has(antenna.id)
            )
        ).map(deepClone),
    };
}

export function createTemplateFromScenarioItem(
    scenarioData: ScenarioData,
    itemId: string,
    options: CreateTemplateOptions = {}
): AssetLibraryTemplate | null {
    const waveform = scenarioData.waveforms.find((item) => item.id === itemId);
    if (waveform) {
        return {
            ...createBaseTemplate(waveform, scenarioData, options),
            kind: 'waveform',
            payload: deepClone(waveform),
        };
    }

    const timing = scenarioData.timings.find((item) => item.id === itemId);
    if (timing) {
        return {
            ...createBaseTemplate(timing, scenarioData, options),
            kind: 'timing',
            payload: deepClone(timing),
        };
    }

    const antenna = scenarioData.antennas.find((item) => item.id === itemId);
    if (antenna) {
        return {
            ...createBaseTemplate(antenna, scenarioData, options),
            kind: 'antenna',
            payload: deepClone(antenna),
        };
    }

    const platform = scenarioData.platforms.find((item) => item.id === itemId);
    if (platform) {
        return {
            ...createBaseTemplate(platform, scenarioData, options),
            kind: 'platform',
            payload: sanitizePlatform(platform),
            dependencies: collectPlatformDependencies(platform, scenarioData),
        };
    }

    return null;
}

function cloneWaveform(
    waveform: Waveform,
    scenarioData: ScenarioData,
    baseName = waveform.name
): Waveform {
    return {
        ...deepClone(waveform),
        id: generateSimId('Waveform'),
        name: createUniqueScenarioCopyName(scenarioData, baseName),
    };
}

function cloneTiming(
    timing: Timing,
    scenarioData: ScenarioData,
    baseName = timing.name
): Timing {
    return {
        ...deepClone(timing),
        id: generateSimId('Timing'),
        name: createUniqueScenarioCopyName(scenarioData, baseName),
        noiseEntries: timing.noiseEntries.map((entry) => ({
            ...deepClone(entry),
            id: generateSimId('Timing'),
        })),
    };
}

function cloneAntenna(
    antenna: Antenna,
    scenarioData: ScenarioData,
    baseName = antenna.name
): Antenna {
    return {
        ...deepClone(antenna),
        id: generateSimId('Antenna'),
        name: createUniqueScenarioCopyName(scenarioData, baseName),
    };
}

function clonePlatformWaypointIds(platform: Platform): Platform {
    const cloned = sanitizePlatform(platform);
    cloned.id = generateSimId('Platform');
    cloned.motionPath.waypoints = cloned.motionPath.waypoints.map(
        (waypoint) => ({
            ...waypoint,
            id: generateSimId('Platform'),
        })
    );

    if (cloned.rotation.type === 'path') {
        cloned.rotation.waypoints = cloned.rotation.waypoints.map(
            (waypoint) => ({
                ...waypoint,
                id: generateSimId('Platform'),
            })
        );
    }

    return cloned;
}

function remapReference(
    originalId: string | null,
    idMap: Map<string, string>,
    warnings: string[],
    label: string
): string | null {
    if (!originalId) {
        return null;
    }

    const remapped = idMap.get(originalId);
    if (remapped) {
        return remapped;
    }

    warnings.push(
        `Missing ${label} dependency ${originalId}; reference cleared.`
    );
    return null;
}

function cloneComponent(
    component: PlatformComponent,
    antennaIdMap: Map<string, string>,
    timingIdMap: Map<string, string>,
    waveformIdMap: Map<string, string>,
    warnings: string[],
    existingNames: Set<string>
): PlatformComponent {
    const name = createUniqueName(component.name, existingNames, {
        copy: true,
    });
    existingNames.add(name);

    switch (component.type) {
        case 'monostatic': {
            const txId = generateSimId('Transmitter');
            return {
                ...deepClone(component),
                id: txId,
                txId,
                rxId: generateSimId('Receiver'),
                name,
                antennaId: remapReference(
                    component.antennaId,
                    antennaIdMap,
                    warnings,
                    'antenna'
                ),
                timingId: remapReference(
                    component.timingId,
                    timingIdMap,
                    warnings,
                    'timing'
                ),
                waveformId: remapReference(
                    component.waveformId,
                    waveformIdMap,
                    warnings,
                    'waveform'
                ),
            };
        }
        case 'transmitter':
            return {
                ...deepClone(component),
                id: generateSimId('Transmitter'),
                name,
                antennaId: remapReference(
                    component.antennaId,
                    antennaIdMap,
                    warnings,
                    'antenna'
                ),
                timingId: remapReference(
                    component.timingId,
                    timingIdMap,
                    warnings,
                    'timing'
                ),
                waveformId: remapReference(
                    component.waveformId,
                    waveformIdMap,
                    warnings,
                    'waveform'
                ),
            };
        case 'receiver':
            return {
                ...deepClone(component),
                id: generateSimId('Receiver'),
                name,
                antennaId: remapReference(
                    component.antennaId,
                    antennaIdMap,
                    warnings,
                    'antenna'
                ),
                timingId: remapReference(
                    component.timingId,
                    timingIdMap,
                    warnings,
                    'timing'
                ),
            };
        case 'target':
            return {
                ...deepClone(component),
                id: generateSimId('Target'),
                name,
            };
    }
}

export function cloneTemplateIntoScenarioData(
    scenarioData: ScenarioData,
    template: AssetLibraryTemplate
): { scenarioData: ScenarioData; result: AssetTemplateInsertionResult } {
    const nextScenarioData: ScenarioData = {
        globalParameters: deepClone(scenarioData.globalParameters),
        waveforms: scenarioData.waveforms.map(deepClone),
        timings: scenarioData.timings.map(deepClone),
        antennas: scenarioData.antennas.map(deepClone),
        platforms: scenarioData.platforms.map(sanitizePlatform),
    };
    const warnings: string[] = [];

    switch (template.kind) {
        case 'waveform': {
            const waveform = cloneWaveform(
                template.payload,
                nextScenarioData,
                template.name
            );
            nextScenarioData.waveforms.push(waveform);
            return {
                scenarioData: nextScenarioData,
                result: {
                    insertedItemId: waveform.id,
                    insertedName: waveform.name,
                    warnings,
                },
            };
        }
        case 'timing': {
            const timing = cloneTiming(
                template.payload,
                nextScenarioData,
                template.name
            );
            nextScenarioData.timings.push(timing);
            return {
                scenarioData: nextScenarioData,
                result: {
                    insertedItemId: timing.id,
                    insertedName: timing.name,
                    warnings,
                },
            };
        }
        case 'antenna': {
            const antenna = cloneAntenna(
                template.payload,
                nextScenarioData,
                template.name
            );
            nextScenarioData.antennas.push(antenna);
            return {
                scenarioData: nextScenarioData,
                result: {
                    insertedItemId: antenna.id,
                    insertedName: antenna.name,
                    warnings,
                },
            };
        }
        case 'platform': {
            const dependencies = template.dependencies ?? {
                waveforms: [],
                timings: [],
                antennas: [],
            };
            const waveformIdMap = new Map<string, string>();
            const timingIdMap = new Map<string, string>();
            const antennaIdMap = new Map<string, string>();

            for (const waveformTemplate of dependencies.waveforms) {
                const waveform = cloneWaveform(
                    waveformTemplate,
                    nextScenarioData
                );
                waveformIdMap.set(waveformTemplate.id, waveform.id);
                nextScenarioData.waveforms.push(waveform);
            }

            for (const timingTemplate of dependencies.timings) {
                const timing = cloneTiming(timingTemplate, nextScenarioData);
                timingIdMap.set(timingTemplate.id, timing.id);
                nextScenarioData.timings.push(timing);
            }

            for (const antennaTemplate of dependencies.antennas) {
                const antenna = cloneAntenna(antennaTemplate, nextScenarioData);
                antennaIdMap.set(antennaTemplate.id, antenna.id);
                nextScenarioData.antennas.push(antenna);
            }

            const platform = clonePlatformWaypointIds(template.payload);
            platform.name = createUniqueScenarioCopyName(
                nextScenarioData,
                template.name
            );
            const existingNames = collectScenarioNames(nextScenarioData);
            existingNames.add(platform.name);
            platform.components = template.payload.components.map((component) =>
                cloneComponent(
                    component,
                    antennaIdMap,
                    timingIdMap,
                    waveformIdMap,
                    warnings,
                    existingNames
                )
            );
            nextScenarioData.platforms.push(platform);

            return {
                scenarioData: nextScenarioData,
                result: {
                    insertedItemId: platform.id,
                    insertedName: platform.name,
                    warnings,
                },
            };
        }
    }
}

export function createAssetLibraryFile(
    templates: AssetLibraryTemplate[]
): AssetLibraryFile {
    return {
        schemaVersion: ASSET_LIBRARY_SCHEMA_VERSION,
        templates,
    };
}

export function parseAssetTemplates(data: unknown): AssetLibraryTemplate[] {
    const fileResult = AssetLibraryFileSchema.safeParse(data);
    if (fileResult.success) {
        return fileResult.data.templates;
    }

    const templateResult = AssetLibraryTemplateSchema.safeParse(data);
    if (templateResult.success) {
        return [templateResult.data];
    }

    throw new Error(
        'Asset library JSON does not match the v1 template schema.'
    );
}

export function prepareTemplatesForCatalog(
    templates: AssetLibraryTemplate[],
    timestamp = new Date().toISOString()
): AssetLibraryTemplate[] {
    return templates.map((template) => ({
        ...deepClone(template),
        id: createAssetTemplateId(),
        updatedAt: timestamp,
    }));
}
