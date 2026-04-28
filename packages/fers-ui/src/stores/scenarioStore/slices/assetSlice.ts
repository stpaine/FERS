// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { StateCreator } from 'zustand';
import { defaultAntenna, defaultTiming, defaultWaveform } from '../defaults';
import { generateSimId } from '../idUtils';
import { cleanObject, serializeAntenna, serializeTiming } from '../serializers';
import { enqueueFullSync, enqueueGranularSyncDetached } from '../syncQueue';
import { Antenna, AssetActions, ScenarioStore } from '../types';
import { buildScenarioJson } from './backendSlice';

export const createAssetSlice: StateCreator<
    ScenarioStore,
    [['zustand/immer', never]],
    [],
    AssetActions
> = (set, get) => ({
    addWaveform: () => {
        set((state) => {
            const id = generateSimId('Waveform');
            state.waveforms.push({
                ...defaultWaveform,
                id,
                name: `Waveform ${state.waveforms.length + 1}`,
            });
            state.isDirty = true;
        });
        // libfers has no granular add API for Waveforms — full sync is required.
        void enqueueFullSync(() => buildScenarioJson(get()));
    },
    addTiming: () => {
        set((state) => {
            const id = generateSimId('Timing');
            state.timings.push({
                ...defaultTiming,
                id,
                name: `Timing ${state.timings.length + 1}`,
            });
            state.isDirty = true;
        });
        // libfers has no granular add API for Timings — full sync is required.
        void enqueueFullSync(() => buildScenarioJson(get()));
    },
    addAntenna: () => {
        set((state) => {
            const id = generateSimId('Antenna');
            state.antennas.push({
                ...defaultAntenna,
                id,
                name: `Antenna ${state.antennas.length + 1}`,
            });
            state.isDirty = true;
        });
        // libfers has no granular add API for Antennas — full sync is required.
        void enqueueFullSync(() => buildScenarioJson(get()));
    },
    addNoiseEntry: (timingId) => {
        let touched = false;
        set((state) => {
            const timing = state.timings.find((t) => t.id === timingId);
            if (timing) {
                timing.noiseEntries.push({
                    id: generateSimId('Timing'),
                    alpha: 0,
                    weight: 0,
                });
                state.isDirty = true;
                touched = true;
            }
        });
        if (touched) {
            const timing = get().timings.find((t) => t.id === timingId);
            if (timing) {
                enqueueGranularSyncDetached(
                    'Timing',
                    timing.id,
                    JSON.stringify(cleanObject(serializeTiming(timing)))
                );
            }
        }
    },
    removeNoiseEntry: (timingId, entryId) => {
        let touched = false;
        set((state) => {
            const timing = state.timings.find((t) => t.id === timingId);
            if (timing) {
                const index = timing.noiseEntries.findIndex(
                    (e) => e.id === entryId
                );
                if (index > -1) {
                    timing.noiseEntries.splice(index, 1);
                    state.isDirty = true;
                    touched = true;
                }
            }
        });
        if (touched) {
            const timing = get().timings.find((t) => t.id === timingId);
            if (timing) {
                enqueueGranularSyncDetached(
                    'Timing',
                    timing.id,
                    JSON.stringify(cleanObject(serializeTiming(timing)))
                );
            }
        }
    },
    setAntennaPattern: (antennaId, newPattern) => {
        let touched = false;
        set((state) => {
            const index = state.antennas.findIndex((a) => a.id === antennaId);
            if (index === -1) return;
            delete state.antennaPreviewErrors[antennaId];

            const oldAntenna = state.antennas[index];
            const baseAntenna = {
                id: oldAntenna.id,
                type: oldAntenna.type,
                name: oldAntenna.name,
                efficiency: oldAntenna.efficiency,
                meshScale: oldAntenna.meshScale,
                design_frequency: oldAntenna.design_frequency,
            };

            let newAntennaState: Antenna;

            switch (newPattern) {
                case 'isotropic':
                    newAntennaState = {
                        ...baseAntenna,
                        pattern: 'isotropic',
                    };
                    break;
                case 'sinc':
                    newAntennaState = {
                        ...baseAntenna,
                        pattern: 'sinc',
                        alpha: 1.0,
                        beta: 1.0,
                        gamma: 2.0,
                    };
                    break;
                case 'gaussian':
                    newAntennaState = {
                        ...baseAntenna,
                        pattern: 'gaussian',
                        azscale: 1.0,
                        elscale: 1.0,
                    };
                    break;
                case 'squarehorn':
                case 'parabolic':
                    newAntennaState = {
                        ...baseAntenna,
                        pattern: newPattern,
                        diameter: 0.5,
                    };
                    break;
                case 'xml':
                case 'file':
                    newAntennaState = {
                        ...baseAntenna,
                        pattern: newPattern,
                        filename: '',
                    };
                    break;
            }
            state.antennas[index] = newAntennaState;
            state.isDirty = true;
            touched = true;
        });
        if (touched) {
            const antenna = get().antennas.find((a) => a.id === antennaId);
            if (antenna) {
                enqueueGranularSyncDetached(
                    'Antenna',
                    antenna.id,
                    JSON.stringify(cleanObject(serializeAntenna(antenna)))
                );
            }
        }
    },
});
