// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { StateCreator } from 'zustand';
import { defaultGlobalParameters } from '../defaults';
import { buildHydratedScenarioState, parseScenarioData } from '../hydration';
import {
    cleanObject,
    serializeAntenna,
    serializeComponentInner,
    serializeGlobalParameters,
    serializePlatform,
    serializeTiming,
    serializeWaveform,
} from '../serializers';
import { enqueueFullSync, enqueueGranularSync } from '../syncQueue';
import {
    Antenna,
    Platform,
    ScenarioActions,
    ScenarioStore,
    Timing,
    Waveform,
} from '../types';
import { setPropertyByPath } from '../utils';
import { buildScenarioJson } from './backendSlice';

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

                setPropertyByPath(item, propertyPath, value);
                state.isDirty = true;
                if (item.type === 'Antenna') {
                    delete state.antennaPreviewErrors[item.id];
                }

                if (
                    item.type !== 'Platform' &&
                    item.type !== 'Antenna' &&
                    item.type !== 'Waveform' &&
                    item.type !== 'Timing'
                ) {
                    return;
                }

                const componentMatch = propertyPath.match(
                    /^components\.(\d+)(?:\.|$)/
                );
                if (item.type === 'Platform' && componentMatch) {
                    const compIndex = parseInt(componentMatch[1], 10);
                    const component = (item as Platform).components[compIndex];
                    if (component) {
                        jsonPayload = JSON.stringify(
                            cleanObject(serializeComponentInner(component))
                        );
                        // Map frontend component type to backend expected type string
                        targetItemType =
                            component.type === 'monostatic'
                                ? 'Monostatic'
                                : component.type.charAt(0).toUpperCase() +
                                  component.type.slice(1);
                        targetItemId = component.id;
                    }
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
                } else if (item.type === 'Waveform') {
                    jsonPayload = JSON.stringify(
                        serializeWaveform(item as Waveform)
                    );
                } else if (item.type === 'Timing') {
                    jsonPayload = JSON.stringify(
                        serializeTiming(item as Timing)
                    );
                }
                return;
            }
        });

        if (jsonPayload) {
            void enqueueGranularSync(targetItemType, targetItemId, jsonPayload);
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
            void enqueueFullSync(() => buildScenarioJson(get()));
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
        void enqueueFullSync(() => buildScenarioJson(get()));
    },
    loadScenario: (backendData: unknown) => {
        const scenarioData = parseScenarioData(backendData);
        if (!scenarioData) {
            return;
        }

        set(buildHydratedScenarioState(get(), scenarioData, { isDirty: true }));
    },
});
