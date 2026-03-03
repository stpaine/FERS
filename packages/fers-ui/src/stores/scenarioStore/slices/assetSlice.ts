// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { StateCreator } from 'zustand';
import { ScenarioStore, AssetActions, Antenna } from '../types';
import { defaultWaveform, defaultTiming, defaultAntenna } from '../defaults';
import { generateSimId } from '../idUtils';

export const createAssetSlice: StateCreator<
    ScenarioStore,
    [['zustand/immer', never]],
    [],
    AssetActions
> = (set) => ({
    addWaveform: () =>
        set((state) => {
            const id = generateSimId('Waveform');
            state.waveforms.push({
                ...defaultWaveform,
                id,
                name: `Waveform ${state.waveforms.length + 1}`,
            });
            state.isDirty = true;
        }),
    addTiming: () =>
        set((state) => {
            const id = generateSimId('Timing');
            state.timings.push({
                ...defaultTiming,
                id,
                name: `Timing ${state.timings.length + 1}`,
            });
            state.isDirty = true;
        }),
    addAntenna: () =>
        set((state) => {
            const id = generateSimId('Antenna');
            state.antennas.push({
                ...defaultAntenna,
                id,
                name: `Antenna ${state.antennas.length + 1}`,
            });
            state.isDirty = true;
        }),
    addNoiseEntry: (timingId) =>
        set((state) => {
            const timing = state.timings.find((t) => t.id === timingId);
            if (timing) {
                timing.noiseEntries.push({
                    id: generateSimId('Timing'),
                    alpha: 0,
                    weight: 0,
                });
                state.isDirty = true;
            }
        }),
    removeNoiseEntry: (timingId, entryId) =>
        set((state) => {
            const timing = state.timings.find((t) => t.id === timingId);
            if (timing) {
                const index = timing.noiseEntries.findIndex(
                    (e) => e.id === entryId
                );
                if (index > -1) {
                    timing.noiseEntries.splice(index, 1);
                    state.isDirty = true;
                }
            }
        }),
    setAntennaPattern: (antennaId, newPattern) =>
        set((state) => {
            const index = state.antennas.findIndex((a) => a.id === antennaId);
            if (index === -1) return;

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
        }),
});
