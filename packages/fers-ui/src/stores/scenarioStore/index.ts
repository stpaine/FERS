// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { create } from 'zustand';
import { immer } from 'zustand/middleware/immer';
import { defaultGlobalParameters } from './defaults';
import { createAssetSlice } from './slices/assetSlice';
import { createBackendSlice } from './slices/backendSlice';
import { createPlatformSlice } from './slices/platformSlice';
import { createScenarioSlice } from './slices/scenarioSlice';
import { registerGranularSyncFailureHandler } from './syncQueue';
import { ScenarioStore } from './types';

export * from './types';
export * from './utils';

export const useScenarioStore = create<ScenarioStore>()(
    immer((set, get, store) => ({
        // Initial State
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
        isSimulating: false,
        isGeneratingKml: false,
        isBackendSyncing: false,
        backendVersion: 0,
        scenarioFilePath: null,
        outputDirectory: null,
        antennaPreviewErrors: {},
        errorSnackbar: {
            open: false,
            message: '',
        },
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

        // Slices
        ...createAssetSlice(set, get, store),
        ...createBackendSlice(set, get, store),
        ...createPlatformSlice(set, get, store),
        ...createScenarioSlice(set, get, store),

        // Playback Actions
        togglePlayPause: () =>
            set((state) => ({ isPlaying: !state.isPlaying })),
        setCurrentTime: (time) => {
            const { start, end } = get().globalParameters;
            const newTime =
                typeof time === 'function' ? time(get().currentTime) : time;
            // Clamp time to simulation bounds
            const clampedTime = Math.max(start, Math.min(end, newTime));
            set({ currentTime: clampedTime });
        },
        setTargetPlaybackDuration: (duration) =>
            set({
                targetPlaybackDuration:
                    duration !== null && duration > 0 ? duration : null,
            }),
        setIsSimulating: (isSimulating) => set({ isSimulating }),
        setIsGeneratingKml: (isGeneratingKml) => set({ isGeneratingKml }),

        frameScene: () =>
            set({
                viewControlAction: { type: 'frame', timestamp: Date.now() },
            }),
        focusOnItem: (itemId) =>
            set({
                viewControlAction: {
                    type: 'focus',
                    targetId: itemId,
                    timestamp: Date.now(),
                },
            }),
        toggleFollowItem: (itemId) => {
            const currentAction = get().viewControlAction;
            if (
                currentAction.type === 'follow' &&
                currentAction.targetId === itemId
            ) {
                set({
                    viewControlAction: { type: null, timestamp: Date.now() },
                });
            } else {
                set({
                    viewControlAction: {
                        type: 'follow',
                        targetId: itemId,
                        timestamp: Date.now(),
                    },
                });
            }
        },
        clearViewControlAction: () =>
            set((state) => {
                if (state.viewControlAction.type !== 'follow') {
                    return {
                        viewControlAction: {
                            type: null,
                            timestamp: state.viewControlAction.timestamp,
                        },
                    };
                }
                return {};
            }),
        toggleLayer: (layer) =>
            set((state) => {
                state.visibility[layer] = !state.visibility[layer];
            }),

        // Error Actions
        showError: (message) => set({ errorSnackbar: { open: true, message } }),
        hideError: () =>
            set((state) => ({
                errorSnackbar: { ...state.errorSnackbar, open: false },
            })),
        setAntennaPreviewError: (antennaId, message) =>
            set((state) => {
                state.antennaPreviewErrors[antennaId] = message;
            }),
        clearAntennaPreviewError: (antennaId) =>
            set((state) => {
                delete state.antennaPreviewErrors[antennaId];
            }),
    }))
);

function formatSyncError(error: unknown): string {
    return error instanceof Error ? error.message : String(error);
}

registerGranularSyncFailureHandler(async ({ itemType, itemId, error }) => {
    const syncMessage = `Sync error for ${itemType} ${itemId}: ${formatSyncError(error)}. Reverting to backend state.`;
    useScenarioStore.getState().showError(syncMessage);

    try {
        await useScenarioStore.getState().fetchFromBackend();
    } catch (reloadError) {
        useScenarioStore
            .getState()
            .showError(
                `${syncMessage} Failed to reload backend state: ${formatSyncError(reloadError)}`
            );
    }
});
