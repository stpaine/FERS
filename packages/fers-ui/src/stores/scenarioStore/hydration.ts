// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { ScenarioDataSchema } from '../scenarioSchema';
import { createDefaultPlatform } from './defaults';
import {
    generateSimId,
    normalizeSimId,
    reserveSimId,
    seedSimIdCounters,
} from './idUtils';
import {
    Antenna,
    FixedRotation,
    GlobalParameters,
    MotionPath,
    Platform,
    PlatformComponent,
    RotationPath,
    ScenarioData,
    ScenarioState,
    Timing,
    Waveform,
} from './types';

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
    fmcw_mode?: object;
    schedule?: BackendSchedulePeriod[];
    rcs?: { type: 'isotropic' | 'file'; value?: number; filename?: string };
    model?: { type: 'constant' | 'chisquare' | 'gamma'; k?: number };
}

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
    fmcw_linear_chirp?: {
        direction: 'up' | 'down';
        chirp_bandwidth: number;
        chirp_duration: number;
        chirp_period: number;
        start_frequency_offset?: number | null;
        chirp_count?: number | null;
    };
    fmcw_triangle?: {
        chirp_bandwidth: number;
        chirp_duration: number;
        start_frequency_offset?: number | null;
        triangle_count?: number | null;
    };
}

type HydratedScenarioState = Pick<
    ScenarioState,
    | 'globalParameters'
    | 'waveforms'
    | 'timings'
    | 'antennas'
    | 'platforms'
    | 'selectedItemId'
    | 'selectedComponentId'
    | 'isDirty'
    | 'antennaPreviewErrors'
    | 'currentTime'
>;

type HydrationOptions = {
    isDirty: boolean;
    preserveSelection?: boolean;
    preserveCurrentTime?: boolean;
};

function hasScenarioItem(scenarioData: ScenarioData, itemId: string): boolean {
    if (itemId === 'global-parameters') {
        return true;
    }

    return [
        ...scenarioData.waveforms,
        ...scenarioData.timings,
        ...scenarioData.antennas,
        ...scenarioData.platforms,
    ].some((item) => item.id === itemId);
}

function resolveSelection(
    currentState: ScenarioState,
    scenarioData: ScenarioData,
    preserveSelection: boolean
): Pick<HydratedScenarioState, 'selectedItemId' | 'selectedComponentId'> {
    if (!preserveSelection) {
        return {
            selectedItemId: null,
            selectedComponentId: null,
        };
    }

    const { selectedItemId, selectedComponentId } = currentState;
    if (selectedComponentId) {
        for (const platform of scenarioData.platforms) {
            const component = platform.components.find(
                (candidate) => candidate.id === selectedComponentId
            );
            if (component) {
                return {
                    selectedItemId: platform.id,
                    selectedComponentId: component.id,
                };
            }
        }
    }

    if (selectedItemId && hasScenarioItem(scenarioData, selectedItemId)) {
        return {
            selectedItemId,
            selectedComponentId: null,
        };
    }

    return {
        selectedItemId: null,
        selectedComponentId: null,
    };
}

function clampCurrentTime(
    currentTime: number,
    globalParameters: GlobalParameters
): number {
    return Math.max(
        globalParameters.start,
        Math.min(globalParameters.end, currentTime)
    );
}

export function buildHydratedScenarioState(
    currentState: ScenarioState,
    scenarioData: ScenarioData,
    options: HydrationOptions
): HydratedScenarioState {
    const selection = resolveSelection(
        currentState,
        scenarioData,
        options.preserveSelection ?? false
    );

    return {
        ...scenarioData,
        ...selection,
        isDirty: options.isDirty,
        antennaPreviewErrors: {},
        currentTime: options.preserveCurrentTime
            ? clampCurrentTime(
                  currentState.currentTime,
                  scenarioData.globalParameters
              )
            : scenarioData.globalParameters.start,
    };
}

export function parseScenarioData(backendData: unknown): ScenarioData | null {
    try {
        if (typeof backendData !== 'object' || backendData === null) {
            console.error('Invalid backend data format: not an object.');
            return null;
        }

        const data =
            'simulation' in backendData &&
            typeof (backendData as { simulation: unknown }).simulation ===
                'object' &&
            (backendData as { simulation: unknown }).simulation !== null
                ? ((backendData as { simulation: object }).simulation as Record<
                      string,
                      unknown
                  >)
                : (backendData as Record<string, unknown>);

        const nameToIdMap = new Map<string, string>();
        const normalizeIdOrNull = (value: unknown) => normalizeSimId(value);
        const normalizeRefId = (value: unknown) => {
            const id = normalizeIdOrNull(value);
            if (!id || id === '0') return null;
            return id;
        };

        const params = (data.parameters as Record<string, unknown>) || {};
        const globalParameters: GlobalParameters = {
            id: 'global-parameters',
            type: 'GlobalParameters',
            rotationAngleUnit:
                (params.rotationangleunit as 'deg' | 'rad') ?? 'deg',
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
                zone: (params.coordinatesystem as Record<string, number>)?.zone,
                hemisphere: (
                    params.coordinatesystem as Record<string, 'N' | 'S'>
                )?.hemisphere,
            },
        };

        const waveforms: Waveform[] = (
            (data.waveforms as BackendWaveform[]) || []
        ).map((w) => {
            const waveformId =
                normalizeIdOrNull(w.id) ?? generateSimId('Waveform');
            const commonWaveform = {
                id: waveformId,
                type: 'Waveform' as const,
                name: w.name,
                power: w.power,
                carrier_frequency: w.carrier_frequency,
            };
            const waveform: Waveform = w.fmcw_linear_chirp
                ? {
                      ...commonWaveform,
                      waveformType: 'fmcw_linear_chirp',
                      direction: w.fmcw_linear_chirp.direction,
                      chirp_bandwidth: w.fmcw_linear_chirp.chirp_bandwidth,
                      chirp_duration: w.fmcw_linear_chirp.chirp_duration,
                      chirp_period: w.fmcw_linear_chirp.chirp_period,
                      start_frequency_offset:
                          w.fmcw_linear_chirp.start_frequency_offset ?? null,
                      chirp_count: w.fmcw_linear_chirp.chirp_count ?? null,
                  }
                : w.fmcw_triangle
                  ? {
                        ...commonWaveform,
                        waveformType: 'fmcw_triangle',
                        chirp_bandwidth: w.fmcw_triangle.chirp_bandwidth,
                        chirp_duration: w.fmcw_triangle.chirp_duration,
                        start_frequency_offset:
                            w.fmcw_triangle.start_frequency_offset ?? null,
                        triangle_count: w.fmcw_triangle.triangle_count ?? null,
                    }
                  : w.cw
                    ? {
                          ...commonWaveform,
                          waveformType: 'cw',
                      }
                    : w.pulsed_from_file
                      ? {
                            ...commonWaveform,
                            waveformType: 'pulsed_from_file',
                            filename: w.pulsed_from_file.filename ?? '',
                        }
                      : (() => {
                            throw new Error(
                                `Unsupported waveform type for '${w.name}'.`
                            );
                        })();
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
                          : cData.fmcw_mode
                            ? 'fmcw'
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
                                window_length: pulsed?.window_length ?? null,
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
                                window_length: pulsed?.window_length ?? null,
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
                    c.type === 'monostatic' ? [c.id, c.txId, c.rxId] : [c.id]
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

        const result = ScenarioDataSchema.safeParse(transformedScenario);
        if (!result.success) {
            console.error(
                'Data validation failed after loading from backend:',
                result.error.flatten()
            );
            return null;
        }

        return result.data;
    } catch (error) {
        console.error(
            'An unexpected error occurred while loading the scenario:',
            error
        );
        return null;
    }
}
