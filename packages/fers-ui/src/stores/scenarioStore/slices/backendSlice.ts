// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { invoke } from '@tauri-apps/api/core';
import { StateCreator } from 'zustand';
import { buildHydratedScenarioState, parseScenarioData } from '../hydration';
import {
    serializeAntenna,
    serializeGlobalParameters,
    serializePlatform,
    serializeTiming,
    serializeWaveform,
} from '../serializers';
import { enqueueFullSync } from '../syncQueue';
import { BackendActions, ScenarioState, ScenarioStore } from '../types';

/**
 * Build the full scenario JSON payload expected by the `update_scenario_from_json`
 * Tauri command. Extracted from `syncBackend` so the sync queue can capture a
 * snapshot at task-execution time rather than enqueue time.
 */
export function buildScenarioJson(state: ScenarioState): string {
    const { globalParameters, waveforms, timings, antennas, platforms } = state;
    const scenarioJson = {
        simulation: {
            name: globalParameters.simulation_name,
            parameters: serializeGlobalParameters(globalParameters),
            waveforms: waveforms.map(serializeWaveform),
            timings: timings.map(serializeTiming),
            antennas: antennas.map(serializeAntenna),
            platforms: platforms.map(serializePlatform),
        },
    };
    return JSON.stringify(scenarioJson, null, 2);
}

export const createBackendSlice: StateCreator<
    ScenarioStore,
    [['zustand/immer', never]],
    [],
    BackendActions
> = (set, get) => ({
    syncBackend: async () => {
        set({ isBackendSyncing: true });
        try {
            await enqueueFullSync(() => buildScenarioJson(get()));
            set((state) => {
                state.isBackendSyncing = false;
                state.backendVersion += 1;
            });
        } catch (error) {
            set({ isBackendSyncing: false });
            throw error;
        }
    },
    fetchFromBackend: async () => {
        try {
            const jsonState = await invoke<string>('get_scenario_as_json');
            const parsedJson = JSON.parse(jsonState);
            const scenarioData = parseScenarioData(parsedJson);
            if (!scenarioData) {
                throw new Error('Failed to hydrate scenario from backend JSON');
            }

            set(
                buildHydratedScenarioState(get(), scenarioData, {
                    isDirty: false,
                    preserveSelection: true,
                    preserveCurrentTime: true,
                })
            );
        } catch (error) {
            console.error('Failed to fetch state from backend:', error);
            throw error;
        }
    },
});
