// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { StateCreator } from 'zustand';
import { cloneTemplateIntoScenarioData } from '../../assetTemplates';
import { useSimulationProgressStore } from '../../simulationProgressStore';
import { defaultGlobalParameters } from '../defaults';
import { buildHydratedScenarioState, parseScenarioData } from '../hydration';
import {
    createUniqueScenarioName,
    getComponentIdentityIds,
} from '../nameUtils';
import {
    serializeAntenna,
    serializeGlobalParameters,
    serializePlatform,
    serializeTiming,
} from '../serializers';
import {
    enqueueFullSyncDetached,
    enqueueGranularSyncDetached,
} from '../syncQueue';
import {
    Antenna,
    GlobalParameters,
    Platform,
    ScenarioActions,
    ScenarioStore,
    Timing,
} from '../types';
import { setPropertyByPath } from '../utils';
import { buildScenarioJson } from './backendSlice';

function convertRotationValue(
    value: number,
    fromUnit: GlobalParameters['rotationAngleUnit'],
    toUnit: GlobalParameters['rotationAngleUnit']
): number {
    if (fromUnit === toUnit) {
        return value;
    }
    return fromUnit === 'deg'
        ? (value * Math.PI) / 180
        : (value * 180) / Math.PI;
}

export const createScenarioSlice: StateCreator<
    ScenarioStore,
    [['zustand/immer', never]],
    [],
    ScenarioActions
> = (set, get) => ({
    setScenarioFilePath: (path) => set({ scenarioFilePath: path }),
    setOutputDirectory: (dir) => set({ outputDirectory: dir }),

    selectItem: (itemId) =>
        set((state) => {
            if (!itemId) {
                state.selectedItemId = null;
                state.selectedComponentId = null;
                return;
            }

            if (itemId === 'global-parameters') {
                state.selectedItemId = itemId;
                state.selectedComponentId = null;
                return;
            }

            const collections = [
                'waveforms',
                'timings',
                'antennas',
                'platforms',
            ] as const;

            for (const key of collections) {
                if (state[key].some((i: { id: string }) => i.id === itemId)) {
                    state.selectedItemId = itemId;
                    state.selectedComponentId = null;
                    return;
                }
            }

            for (const platform of state.platforms) {
                const component = platform.components.find(
                    (c) => c.id === itemId
                );
                if (component) {
                    state.selectedItemId = platform.id;
                    state.selectedComponentId = component.id;
                    return;
                }
                const monostatic = platform.components.find(
                    (c) =>
                        c.type === 'monostatic' &&
                        (c.txId === itemId || c.rxId === itemId)
                );
                if (monostatic) {
                    state.selectedItemId = platform.id;
                    state.selectedComponentId = monostatic.id;
                    return;
                }
            }

            state.selectedItemId = itemId;
            state.selectedComponentId = null;
        }),
    updateItem: (itemId, propertyPath, value) => {
        let targetItemType = '';
        let targetItemId = '';
        let jsonPayload: string | null = null;
        let requiresFullSync = false;

        set((state) => {
            if (itemId === 'global-parameters') {
                setPropertyByPath(state.globalParameters, propertyPath, value);

                // Ensure currentTime remains within valid bounds if start/end parameters are updated
                if (propertyPath === 'start' && typeof value === 'number') {
                    state.currentTime = Math.max(state.currentTime, value);
                } else if (
                    propertyPath === 'end' &&
                    typeof value === 'number'
                ) {
                    state.currentTime = Math.min(state.currentTime, value);
                }

                state.isDirty = true;

                if (
                    ['start', 'end', 'rate', 'oversample_ratio'].includes(
                        propertyPath
                    )
                ) {
                    requiresFullSync = true;
                    return;
                }

                targetItemType = 'GlobalParameters';
                targetItemId = itemId;
                jsonPayload = JSON.stringify(
                    serializeGlobalParameters(state.globalParameters)
                );
                return;
            }

            const collections = [
                'waveforms',
                'timings',
                'antennas',
                'platforms',
            ] as const;

            for (const key of collections) {
                const item = state[key].find(
                    (i: { id: string }) => i.id === itemId
                );
                if (!item) continue;

                let nextValue = value;
                if (typeof value === 'string') {
                    if (propertyPath === 'name') {
                        nextValue = createUniqueScenarioName(state, value, [
                            item.id,
                        ]);
                    } else if (item.type === 'Platform') {
                        const componentNameMatch = propertyPath.match(
                            /^components\.(\d+)\.name$/
                        );
                        if (componentNameMatch) {
                            const component =
                                item.components[Number(componentNameMatch[1])];
                            if (component) {
                                nextValue = createUniqueScenarioName(
                                    state,
                                    value,
                                    getComponentIdentityIds(component)
                                );
                            }
                        }
                    }
                }

                setPropertyByPath(item, propertyPath, nextValue);
                state.isDirty = true;
                if (item.type === 'Antenna') {
                    delete state.antennaPreviewErrors[item.id];
                }

                if (item.type === 'Waveform') {
                    requiresFullSync = true;
                    return;
                }

                if (
                    item.type !== 'Platform' &&
                    item.type !== 'Antenna' &&
                    item.type !== 'Timing'
                ) {
                    return;
                }

                const componentMatch = propertyPath.match(
                    /^components\.(\d+)(?:\.|$)/
                );
                if (item.type === 'Platform' && componentMatch) {
                    // Backend full-scenario parsing skips incomplete child
                    // components until required references or file paths exist.
                    // Rebuild from the snapshot so draft components become real
                    // as soon as their authoring state is complete.
                    requiresFullSync = true;
                    return;
                }

                targetItemType = item.type;
                targetItemId = item.id;
                if (item.type === 'Platform') {
                    jsonPayload = JSON.stringify(
                        serializePlatform(item as Platform)
                    );
                } else if (item.type === 'Antenna') {
                    jsonPayload = JSON.stringify(
                        serializeAntenna(item as Antenna)
                    );
                } else if (item.type === 'Timing') {
                    jsonPayload = JSON.stringify(
                        serializeTiming(item as Timing)
                    );
                }
                return;
            }
        });

        if (requiresFullSync) {
            enqueueFullSyncDetached(() => buildScenarioJson(get()));
        } else if (jsonPayload) {
            enqueueGranularSyncDetached(
                targetItemType,
                targetItemId,
                jsonPayload
            );
        }
    },
    setRotationAngleUnit: (unit, convertExisting) => {
        const previousUnit = get().globalParameters.rotationAngleUnit;
        if (previousUnit === unit) {
            return;
        }

        set((state) => {
            state.globalParameters.rotationAngleUnit = unit;

            if (convertExisting) {
                for (const platform of state.platforms) {
                    if (platform.rotation.type === 'fixed') {
                        platform.rotation.startAzimuth = convertRotationValue(
                            platform.rotation.startAzimuth,
                            previousUnit,
                            unit
                        );
                        platform.rotation.startElevation = convertRotationValue(
                            platform.rotation.startElevation,
                            previousUnit,
                            unit
                        );
                        platform.rotation.azimuthRate = convertRotationValue(
                            platform.rotation.azimuthRate,
                            previousUnit,
                            unit
                        );
                        platform.rotation.elevationRate = convertRotationValue(
                            platform.rotation.elevationRate,
                            previousUnit,
                            unit
                        );
                    } else {
                        for (const waypoint of platform.rotation.waypoints) {
                            waypoint.azimuth = convertRotationValue(
                                waypoint.azimuth,
                                previousUnit,
                                unit
                            );
                            waypoint.elevation = convertRotationValue(
                                waypoint.elevation,
                                previousUnit,
                                unit
                            );
                        }
                    }
                }
            }

            for (const platform of state.platforms) {
                delete platform.rotationPathPoints;
            }
            state.isDirty = true;
        });

        enqueueFullSyncDetached(() => buildScenarioJson(get()));
        for (const platform of get().platforms) {
            void get().fetchPlatformPath(platform.id);
        }
    },
    removeItem: (itemId) => {
        let removed = false;
        set((state) => {
            if (!itemId) return;
            const collections = [
                'waveforms',
                'timings',
                'antennas',
                'platforms',
            ] as const;
            for (const key of collections) {
                const index = state[key].findIndex(
                    (item: { id: string }) => item.id === itemId
                );
                if (index > -1) {
                    if (key === 'antennas') {
                        delete state.antennaPreviewErrors[itemId];
                    }
                    state[key].splice(index, 1);
                    if (state.selectedItemId === itemId) {
                        state.selectedItemId = null;
                    }
                    if (key === 'platforms') {
                        state.selectedComponentId = null;
                    }
                    state.isDirty = true;
                    removed = true;
                    return;
                }
            }
        });
        if (removed) {
            // libfers has no granular remove API — full sync is required.
            enqueueFullSyncDetached(() => buildScenarioJson(get()));
        }
    },
    resetScenario: () => {
        set({
            globalParameters: defaultGlobalParameters,
            waveforms: [],
            timings: [],
            antennas: [],
            platforms: [],
            selectedItemId: null,
            selectedComponentId: null,
            isDirty: false,
            antennaPreviewErrors: {},
            currentTime: defaultGlobalParameters.start,
        });
        useSimulationProgressStore.getState().clearSimulationProgress();
        useSimulationProgressStore.setState({ isSimulating: false });
        enqueueFullSyncDetached(() => buildScenarioJson(get()));
    },
    loadScenario: (backendData: unknown) => {
        const scenarioData = parseScenarioData(backendData);
        if (!scenarioData) {
            return;
        }

        set(buildHydratedScenarioState(get(), scenarioData, { isDirty: true }));
    },
    insertAssetTemplate: (template) => {
        const { scenarioData, result } = cloneTemplateIntoScenarioData(
            get(),
            template
        );

        set((state) => {
            state.waveforms = scenarioData.waveforms;
            state.timings = scenarioData.timings;
            state.antennas = scenarioData.antennas;
            state.platforms = scenarioData.platforms;
            state.selectedItemId = result.insertedItemId;
            state.selectedComponentId = null;
            state.isDirty = true;
        });

        // v1 loads library templates by cloning them into the active scenario.
        // Replacing selected items or prompting for placement each time remain
        // future options; nested component-only templates are intentionally out
        // of scope until whole-platform reuse proves insufficient.
        result.warnings.forEach((warning) => get().showWarning(warning));
        enqueueFullSyncDetached(() => buildScenarioJson(get()));

        return result;
    },
});
