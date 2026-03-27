// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { StateCreator } from 'zustand';
import { invoke } from '@tauri-apps/api/core';
import {
    ScenarioStore,
    ScenarioActions,
    GlobalParameters,
    Waveform,
    Timing,
    Antenna,
    Platform,
    MotionPath,
    FixedRotation,
    RotationPath,
    PlatformComponent,
    ScenarioData,
} from '../types';
import { createDefaultPlatform, defaultGlobalParameters } from '../defaults';
import { setPropertyByPath } from '../utils';
import { ScenarioDataSchema } from '../../scenarioSchema';
import {
    generateSimId,
    normalizeSimId,
    seedSimIdCounters,
    reserveSimId,
} from '../idUtils';

// Define interfaces for the expected backend data structure.
interface BackendObjectWithName {
    name: string;
    [key: string]: unknown;
}

interface BackendPulsedMode {
    prf?: number;
    window_skip?: number;
    window_length?: number;
}

interface BackendSchedulePeriod {
    start: number;
    end: number;
}

interface BackendPlatformComponentData {
    name: string;
    id?: string | number;
    tx_id?: string | number;
    rx_id?: string | number;
    antenna?: string | number;
    timing?: string | number;
    waveform?: string | number;
    noise_temp?: number | null;
    nodirect?: boolean;
    nopropagationloss?: boolean;
    pulsed_mode?: BackendPulsedMode;
    cw_mode?: object;
    schedule?: BackendSchedulePeriod[];
    rcs?: { type: 'isotropic' | 'file'; value?: number; filename?: string };
    model?: { type: 'constant' | 'chisquare' | 'gamma'; k?: number };
}

// Backend waypoint types (frontend type minus 'id')
interface BackendPositionWaypoint {
    x: number;
    y: number;
    altitude: number;
    time: number;
}

interface BackendRotationWaypoint {
    azimuth: number;
    elevation: number;
    time: number;
}

interface BackendPlatform {
    name: string;
    id?: string | number;
    motionpath?: {
        interpolation: 'static' | 'linear' | 'cubic';
        positionwaypoints?: BackendPositionWaypoint[];
    };
    fixedrotation?: {
        startazimuth: number;
        startelevation: number;
        azimuthrate: number;
        elevationrate: number;
    };
    rotationpath?: {
        interpolation: 'static' | 'linear' | 'cubic';
        rotationwaypoints?: BackendRotationWaypoint[];
    };
    components?: Record<string, BackendPlatformComponentData>[];
}

interface BackendWaveform {
    name: string;
    id?: string | number;
    power: number;
    carrier_frequency: number;
    cw?: object;
    pulsed_from_file?: {
        filename: string;
    };
}

export const createScenarioSlice: StateCreator<
    ScenarioStore,
    [['zustand/immer', never]],
    [],
    ScenarioActions
> = (set) => ({
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
    updateItem: (itemId, propertyPath, value) =>
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
                if (item) {
                    setPropertyByPath(item, propertyPath, value);
                    state.isDirty = true;

                    // Immediately trigger granular sync for high-performance features
                    // except for complex topology changes which rely on full sync.
                    if (
                        item.type === 'Platform' ||
                        item.type === 'Antenna' ||
                        item.type === 'Waveform'
                    ) {
                        try {
                            let jsonPayload: string | null = null;
                            if (item.type === 'Platform') {
                                const p = item as Platform;
                                const backendRotation: Record<string, unknown> =
                                    {};
                                if (p.rotation.type === 'fixed') {
                                    const r = p.rotation;
                                    backendRotation.fixedrotation = {
                                        interpolation: 'constant',
                                        startazimuth: r.startAzimuth,
                                        startelevation: r.startElevation,
                                        azimuthrate: r.azimuthRate,
                                        elevationrate: r.elevationRate,
                                    };
                                } else {
                                    const r = p.rotation;
                                    backendRotation.rotationpath = {
                                        interpolation: r.interpolation,
                                        rotationwaypoints: r.waypoints.map(
                                            // eslint-disable-next-line @typescript-eslint/no-unused-vars
                                            ({ id, ...wp }) => wp
                                        ),
                                    };
                                }
                                jsonPayload = JSON.stringify({
                                    id: p.id,
                                    name: p.name,
                                    motionpath: {
                                        interpolation:
                                            p.motionPath.interpolation,
                                        positionwaypoints:
                                            p.motionPath.waypoints.map(
                                                // eslint-disable-next-line @typescript-eslint/no-unused-vars
                                                ({ id, ...wp }) => wp
                                            ),
                                    },
                                    ...backendRotation,
                                });
                            } else if (item.type === 'Antenna') {
                                const a = item as Antenna;
                                // eslint-disable-next-line @typescript-eslint/no-unused-vars
                                const { type, meshScale, ...rest } = a;
                                jsonPayload = JSON.stringify(rest);
                            } else if (item.type === 'Waveform') {
                                const w = item as Waveform;
                                const waveformContent =
                                    w.waveformType === 'cw'
                                        ? { cw: {} }
                                        : {
                                              pulsed_from_file: {
                                                  filename: w.filename,
                                              },
                                          };
                                jsonPayload = JSON.stringify({
                                    id: w.id,
                                    name: w.name,
                                    power: w.power,
                                    carrier_frequency: w.carrier_frequency,
                                    ...waveformContent,
                                });
                            }

                            if (jsonPayload) {
                                void invoke('update_item_from_json', {
                                    itemType: item.type,
                                    itemId: item.id,
                                    json: jsonPayload,
                                });
                            }
                        } catch (e) {
                            console.error('Granular sync failed:', e);
                        }
                    }

                    return;
                }
            }
        }),
    removeItem: (itemId) =>
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
                    state[key].splice(index, 1);
                    if (state.selectedItemId === itemId) {
                        state.selectedItemId = null;
                    }
                    if (key === 'platforms') {
                        state.selectedComponentId = null;
                    }
                    state.isDirty = true;
                    return;
                }
            }
        }),
    resetScenario: () =>
        set({
            globalParameters: defaultGlobalParameters,
            waveforms: [],
            timings: [],
            antennas: [],
            platforms: [],
            selectedItemId: null,
            selectedComponentId: null,
            isDirty: false,
            currentTime: defaultGlobalParameters.start,
        }),
    loadScenario: (backendData: unknown) => {
        try {
            if (typeof backendData !== 'object' || backendData === null) {
                console.error('Invalid backend data format: not an object.');
                return;
            }

            const data =
                'simulation' in backendData &&
                typeof (backendData as { simulation: unknown }).simulation ===
                    'object' &&
                (backendData as { simulation: unknown }).simulation !== null
                    ? ((backendData as { simulation: object })
                          .simulation as Record<string, unknown>)
                    : (backendData as Record<string, unknown>);

            const nameToIdMap = new Map<string, string>();
            const normalizeIdOrNull = (value: unknown) => normalizeSimId(value);
            const normalizeRefId = (value: unknown) => {
                const id = normalizeIdOrNull(value);
                if (!id || id === '0') return null;
                return id;
            };

            // 1. Parameters
            const params = (data.parameters as Record<string, unknown>) || {};
            const globalParameters: GlobalParameters = {
                id: 'global-parameters',
                type: 'GlobalParameters',
                simulation_name: (data.name as string) || 'FERS Simulation',
                start: (params.starttime as number) ?? 0.0,
                end: (params.endtime as number) ?? 10.0,
                rate: (params.rate as number) ?? 10000.0,
                simSamplingRate: (params.simSamplingRate as number) ?? null,
                c: (params.c as number) ?? 299792458.0,
                random_seed: (params.randomseed as number) ?? null,
                adc_bits: (params.adc_bits as number) ?? 12,
                oversample_ratio: (params.oversample as number) ?? 1,
                origin: {
                    latitude:
                        ((params.origin as Record<string, number>)
                            ?.latitude as number) ?? -33.957652,
                    longitude:
                        ((params.origin as Record<string, number>)
                            ?.longitude as number) ?? 18.4611991,
                    altitude:
                        ((params.origin as Record<string, number>)
                            ?.altitude as number) ?? 111.01,
                },
                coordinateSystem: {
                    frame:
                        ((params.coordinatesystem as Record<string, string>)
                            ?.frame as GlobalParameters['coordinateSystem']['frame']) ??
                        'ENU',
                    zone: (params.coordinatesystem as Record<string, number>)
                        ?.zone,
                    hemisphere: (
                        params.coordinatesystem as Record<string, 'N' | 'S'>
                    )?.hemisphere,
                },
            };

            // 2. Assets (and build name-to-id map)
            const waveforms: Waveform[] = (
                (data.waveforms as BackendWaveform[]) || []
            ).map((w) => {
                const waveformType = w.cw
                    ? ('cw' as const)
                    : ('pulsed_from_file' as const);
                const filename = w.pulsed_from_file?.filename ?? '';

                const waveformId =
                    normalizeIdOrNull(w.id) ?? generateSimId('Waveform');
                const waveform: Waveform = {
                    id: waveformId,
                    type: 'Waveform',
                    name: w.name,
                    waveformType,
                    power: w.power,
                    carrier_frequency: w.carrier_frequency,
                    filename,
                };
                nameToIdMap.set(waveform.name, waveform.id);
                reserveSimId(waveformId);
                return waveform;
            });

            const timings: Timing[] = (
                (data.timings as BackendObjectWithName[]) || []
            ).map((t) => {
                const timingId =
                    normalizeIdOrNull((t as { id?: unknown }).id) ??
                    generateSimId('Timing');
                const timing = {
                    ...t,
                    id: timingId,
                    type: 'Timing' as const,
                    freqOffset: (t.freq_offset as number) ?? null,
                    randomFreqOffsetStdev:
                        (t.random_freq_offset_stdev as number) ?? null,
                    phaseOffset: (t.phase_offset as number) ?? null,
                    randomPhaseOffsetStdev:
                        (t.random_phase_offset_stdev as number) ?? null,
                    noiseEntries: Array.isArray(t.noise_entries)
                        ? t.noise_entries.map((item) => ({
                              ...(item as object),
                              id: generateSimId('Timing'),
                          }))
                        : [],
                };
                timing.noiseEntries.forEach((entry) => reserveSimId(entry.id));
                nameToIdMap.set(timing.name, timing.id);
                reserveSimId(timingId);
                return timing as Timing;
            });

            const antennas: Antenna[] = (
                (data.antennas as BackendObjectWithName[]) || []
            ).map((a) => {
                const antennaId =
                    normalizeIdOrNull((a as { id?: unknown }).id) ??
                    generateSimId('Antenna');
                const antenna = {
                    ...a,
                    id: antennaId,
                    type: 'Antenna' as const,
                };
                nameToIdMap.set(antenna.name, antenna.id);
                reserveSimId(antennaId);
                return antenna as Antenna;
            });

            // 3. Platforms
            const platforms: Platform[] = (
                (data.platforms as BackendPlatform[]) || []
            ).map((p): Platform => {
                const platformId =
                    normalizeIdOrNull(p.id) ?? generateSimId('Platform');
                reserveSimId(platformId);

                const motionPath: MotionPath = {
                    interpolation: p.motionpath?.interpolation ?? 'static',
                    waypoints: (p.motionpath?.positionwaypoints ?? []).map(
                        (item) => ({
                            ...item,
                            id: generateSimId('Platform'),
                        })
                    ),
                };
                motionPath.waypoints.forEach((wp) => reserveSimId(wp.id));

                let rotation: FixedRotation | RotationPath;
                if (p.fixedrotation) {
                    rotation = {
                        type: 'fixed',
                        startAzimuth: p.fixedrotation.startazimuth,
                        startElevation: p.fixedrotation.startelevation,
                        azimuthRate: p.fixedrotation.azimuthrate,
                        elevationRate: p.fixedrotation.elevationrate,
                    };
                } else if (p.rotationpath) {
                    rotation = {
                        type: 'path',
                        interpolation: p.rotationpath.interpolation ?? 'static',
                        waypoints: (p.rotationpath.rotationwaypoints || []).map(
                            (item) => ({
                                ...item,
                                id: generateSimId('Platform'),
                            })
                        ),
                    };
                    rotation.waypoints.forEach((wp) => reserveSimId(wp.id));
                } else {
                    rotation = createDefaultPlatform().rotation as
                        | FixedRotation
                        | RotationPath;
                }

                const components: PlatformComponent[] = [];

                if (p.components && Array.isArray(p.components)) {
                    p.components.forEach((compWrapper) => {
                        const cType = Object.keys(compWrapper)[0];
                        const cData = compWrapper[cType];
                        const componentId =
                            normalizeIdOrNull(cData.id) ??
                            (cType === 'transmitter'
                                ? generateSimId('Transmitter')
                                : cType === 'receiver'
                                  ? generateSimId('Receiver')
                                  : cType === 'target'
                                    ? generateSimId('Target')
                                    : generateSimId('Transmitter'));
                        const txId =
                            normalizeIdOrNull(cData.tx_id) ??
                            (cType === 'monostatic'
                                ? generateSimId('Transmitter')
                                : null);
                        const rxId =
                            normalizeIdOrNull(cData.rx_id) ??
                            (cType === 'monostatic'
                                ? generateSimId('Receiver')
                                : null);

                        const radarType = cData.pulsed_mode
                            ? 'pulsed'
                            : cData.cw_mode
                              ? 'cw'
                              : 'pulsed';
                        const pulsed = cData.pulsed_mode;

                        const antennaId =
                            normalizeRefId(cData.antenna) ??
                            nameToIdMap.get(String(cData.antenna ?? '')) ??
                            null;
                        const timingId =
                            normalizeRefId(cData.timing) ??
                            nameToIdMap.get(String(cData.timing ?? '')) ??
                            null;
                        const waveformId =
                            normalizeRefId(cData.waveform) ??
                            nameToIdMap.get(String(cData.waveform ?? '')) ??
                            null;

                        const commonRadar = {
                            antennaId,
                            timingId,
                            schedule: cData.schedule ?? [],
                        };
                        const commonReceiver = {
                            noiseTemperature: cData.noise_temp ?? null,
                            noDirectPaths: cData.nodirect ?? false,
                            noPropagationLoss: cData.nopropagationloss ?? false,
                        };

                        let newComp: PlatformComponent | null = null;

                        switch (cType) {
                            case 'monostatic':
                                newComp = {
                                    id: componentId,
                                    type: 'monostatic',
                                    txId: txId ?? componentId,
                                    rxId: rxId ?? generateSimId('Receiver'),
                                    name: cData.name,
                                    radarType,
                                    window_skip: pulsed?.window_skip ?? null,
                                    window_length:
                                        pulsed?.window_length ?? null,
                                    prf: pulsed?.prf ?? null,
                                    waveformId,
                                    ...commonRadar,
                                    ...commonReceiver,
                                };
                                break;
                            case 'transmitter':
                                newComp = {
                                    id: componentId,
                                    type: 'transmitter',
                                    name: cData.name,
                                    radarType,
                                    prf: pulsed?.prf ?? null,
                                    waveformId,
                                    ...commonRadar,
                                };
                                break;
                            case 'receiver':
                                newComp = {
                                    id: componentId,
                                    type: 'receiver',
                                    name: cData.name,
                                    radarType,
                                    window_skip: pulsed?.window_skip ?? null,
                                    window_length:
                                        pulsed?.window_length ?? null,
                                    prf: pulsed?.prf ?? null,
                                    ...commonRadar,
                                    ...commonReceiver,
                                };
                                break;
                            case 'target':
                                newComp = {
                                    id: componentId,
                                    type: 'target',
                                    name: cData.name,
                                    rcs_type: cData.rcs?.type ?? 'isotropic',
                                    rcs_value: cData.rcs?.value,
                                    rcs_filename: cData.rcs?.filename,
                                    rcs_model: cData.model?.type ?? 'constant',
                                    rcs_k: cData.model?.k,
                                };
                                break;
                        }

                        if (newComp) {
                            components.push(newComp);
                        }
                    });
                }

                return {
                    id: platformId,
                    type: 'Platform',
                    name: p.name,
                    motionPath,
                    rotation,
                    components,
                };
            });

            seedSimIdCounters([
                ...waveforms.map((w) => w.id),
                ...timings.map((t) => t.id),
                ...antennas.map((a) => a.id),
                ...platforms.map((p) => p.id),
                ...platforms.flatMap((p) =>
                    p.components.flatMap((c) =>
                        c.type === 'monostatic'
                            ? [c.id, c.txId, c.rxId]
                            : [c.id]
                    )
                ),
            ]);

            const transformedScenario: ScenarioData = {
                globalParameters,
                waveforms,
                timings,
                antennas,
                platforms,
            };

            // --- Validate and Commit ---
            const result = ScenarioDataSchema.safeParse(transformedScenario);

            if (!result.success) {
                console.error(
                    'Data validation failed after loading from backend:',
                    result.error.flatten()
                );
                return;
            }

            // Update state with the validated and parsed data.
            set({
                ...result.data,
                selectedItemId: null,
                selectedComponentId: null,
                isDirty: true,
                currentTime: result.data.globalParameters.start,
            });
        } catch (error) {
            console.error(
                'An unexpected error occurred while loading the scenario:',
                error
            );
        }
    },
});
