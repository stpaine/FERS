// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { StateCreator } from 'zustand';
import { invoke } from '@tauri-apps/api/core';
import { ScenarioStore, BackendActions } from '../types';
import {
    serializePlatform,
    serializeWaveform,
    serializeTiming,
    serializeAntenna,
    serializeGlobalParameters,
} from '../serializers';

export const createBackendSlice: StateCreator<
    ScenarioStore,
    [['zustand/immer', never]],
    [],
    BackendActions
> = (set, get) => ({
    syncBackend: async () => {
        set({ isBackendSyncing: true });
        const { globalParameters, waveforms, timings, antennas, platforms } =
            get();

        const backendPlatforms = platforms.map(serializePlatform);
        const backendWaveforms = waveforms.map(serializeWaveform);
        const backendTimings = timings.map(serializeTiming);
        const backendAntennas = antennas.map(serializeAntenna);
        const gp_params = serializeGlobalParameters(globalParameters);

        const scenarioJson = {
            simulation: {
                name: globalParameters.simulation_name,
                parameters: gp_params,
                waveforms: backendWaveforms,
                timings: backendTimings,
                antennas: backendAntennas,
                platforms: backendPlatforms,
            },
        };

        try {
            const jsonPayload = JSON.stringify(scenarioJson, null, 2);
            await invoke('update_scenario_from_json', {
                json: jsonPayload,
            });
            console.log('Successfully synced state to backend.');

            set((state) => {
                state.isBackendSyncing = false;
                state.backendVersion += 1;
            });
        } catch (error) {
            console.error('Failed to sync state to backend:', error);
            set({ isBackendSyncing: false });
            throw error;
        }
    },
    fetchFromBackend: async () => {
        try {
            const jsonState = await invoke<string>('get_scenario_as_json');
            const scenarioData = JSON.parse(jsonState);
            get().loadScenario(scenarioData);
        } catch (error) {
            console.error('Failed to fetch state from backend:', error);
            throw error;
        }
    },
});
