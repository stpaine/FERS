// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { omit } from '@/utils/typeUtils';
import {
    Antenna,
    GlobalParameters,
    Platform,
    PlatformComponent,
    Timing,
    Waveform,
} from './types';

type DistributiveOmit<T, K extends PropertyKey> = T extends unknown
    ? Omit<T, Extract<keyof T, K>>
    : never;

type SerializedAntenna =
    | DistributiveOmit<Antenna, 'type' | 'meshScale'>
    | {
          id: Antenna['id'];
          name: Antenna['name'];
          pattern: 'isotropic';
          efficiency: Antenna['efficiency'];
      };

type TargetComponent = Extract<PlatformComponent, { type: 'target' }>;

/**
 * Recursively removes null or undefined values from an object or array.
 */
export const cleanObject = <T>(obj: T): T => {
    if (Array.isArray(obj)) {
        return obj.map(cleanObject) as unknown as T;
    }
    if (obj !== null && typeof obj === 'object') {
        const newObj = {} as Record<string, unknown>;
        for (const [key, value] of Object.entries(obj)) {
            if (value !== null && value !== undefined) {
                newObj[key] = cleanObject(value);
            }
        }
        return newObj as T;
    }
    return obj;
};

const serializeReferenceId = (id: string | null | undefined): string | 0 =>
    typeof id === 'string' && id.trim().length > 0 ? id : 0;

const serializeTargetRcs = (component: TargetComponent) => {
    const filename = component.rcs_filename;
    if (
        component.rcs_type === 'file' &&
        typeof filename === 'string' &&
        filename.trim().length > 0
    ) {
        return {
            type: 'file',
            filename,
        };
    }

    return {
        type: 'isotropic',
        value: component.rcs_value ?? 1,
    };
};

/**
 * Serializes the inner properties of a component (unwrapped).
 * This is the exact format the C++ backend expects for granular updates.
 */
export const serializeComponentInner = (component: PlatformComponent) => {
    const mode = (() => {
        if (!('radarType' in component)) {
            return {};
        }

        switch (component.radarType) {
            case 'pulsed':
                return {
                    pulsed_mode: {
                        prf: component.prf,
                        ...(component.type !== 'transmitter' && {
                            window_skip: component.window_skip,
                            window_length: component.window_length,
                        }),
                    },
                };
            case 'cw':
                return { cw_mode: {} };
            case 'fmcw':
                return { fmcw_mode: {} };
        }
    })();

    switch (component.type) {
        case 'monostatic':
            return {
                tx_id: component.txId,
                rx_id: component.rxId,
                name: component.name,
                ...mode,
                antenna: serializeReferenceId(component.antennaId),
                waveform: serializeReferenceId(component.waveformId),
                timing: serializeReferenceId(component.timingId),
                noise_temp: component.noiseTemperature,
                nodirect: component.noDirectPaths,
                nopropagationloss: component.noPropagationLoss,
                schedule: component.schedule,
            };
        case 'transmitter':
            return {
                id: component.id,
                name: component.name,
                ...mode,
                antenna: serializeReferenceId(component.antennaId),
                waveform: serializeReferenceId(component.waveformId),
                timing: serializeReferenceId(component.timingId),
                schedule: component.schedule,
            };
        case 'receiver':
            return {
                id: component.id,
                name: component.name,
                ...mode,
                antenna: serializeReferenceId(component.antennaId),
                timing: serializeReferenceId(component.timingId),
                noise_temp: component.noiseTemperature,
                nodirect: component.noDirectPaths,
                nopropagationloss: component.noPropagationLoss,
                schedule: component.schedule,
            };
        case 'target': {
            const targetObj: Record<string, unknown> = {
                id: component.id,
                name: component.name,
                rcs: serializeTargetRcs(component),
            };
            if (component.rcs_model !== 'constant') {
                targetObj.model = {
                    type: component.rcs_model,
                    k: component.rcs_k,
                };
            }
            return targetObj;
        }
    }
};

/**
 * Serializes a component wrapped in its type key (e.g., { transmitter: { ... } }).
 * This is the format expected by the full scenario JSON array.
 */
export const serializeComponent = (component: PlatformComponent) => {
    return cleanObject({
        [component.type]: serializeComponentInner(component),
    });
};

export const serializePlatform = (p: Platform) => {
    const { components, motionPath, rotation, ...rest } = p;

    const backendComponents = components.map(serializeComponent);

    const backendRotation: Record<string, unknown> = {};
    if (rotation.type === 'fixed') {
        const r = omit(rotation, 'type');
        backendRotation.fixedrotation = {
            interpolation: 'constant',
            startazimuth: r.startAzimuth,
            startelevation: r.startElevation,
            azimuthrate: r.azimuthRate,
            elevationrate: r.elevationRate,
        };
    } else {
        const r = omit(rotation, 'type');
        backendRotation.rotationpath = {
            interpolation: r.interpolation,
            rotationwaypoints: r.waypoints.map((wp) => omit(wp, 'id')),
        };
    }

    return cleanObject({
        ...rest,
        id: p.id,
        motionpath: {
            interpolation: motionPath.interpolation,
            positionwaypoints: motionPath.waypoints.map((wp) => omit(wp, 'id')),
        },
        ...backendRotation,
        components: backendComponents,
    });
};

export const serializeWaveform = (w: Waveform) => {
    const waveformContent = (() => {
        switch (w.waveformType) {
            case 'cw':
                return { cw: {} };
            case 'pulsed_from_file':
                return { pulsed_from_file: { filename: w.filename } };
            case 'fmcw_linear_chirp':
                return {
                    fmcw_linear_chirp: {
                        direction: w.direction,
                        chirp_bandwidth: w.chirp_bandwidth,
                        chirp_duration: w.chirp_duration,
                        chirp_period: w.chirp_period,
                        ...(w.start_frequency_offset
                            ? {
                                  start_frequency_offset:
                                      w.start_frequency_offset,
                              }
                            : {}),
                        ...(w.chirp_count !== null
                            ? { chirp_count: w.chirp_count }
                            : {}),
                    },
                };
            case 'fmcw_triangle':
                return {
                    fmcw_triangle: {
                        chirp_bandwidth: w.chirp_bandwidth,
                        chirp_duration: w.chirp_duration,
                        ...(w.start_frequency_offset
                            ? {
                                  start_frequency_offset:
                                      w.start_frequency_offset,
                              }
                            : {}),
                        ...(w.triangle_count !== null
                            ? { triangle_count: w.triangle_count }
                            : {}),
                    },
                };
        }
    })();

    return cleanObject({
        id: w.id,
        name: w.name,
        power: w.power,
        carrier_frequency: w.carrier_frequency,
        ...waveformContent,
    });
};

export const serializeTiming = (t: Timing) => {
    const rest = omit(t, 'type');
    const timingObj = {
        ...rest,
        synconpulse: false,
        freq_offset: t.freqOffset,
        random_freq_offset_stdev: t.randomFreqOffsetStdev,
        phase_offset: t.phaseOffset,
        random_phase_offset_stdev: t.randomPhaseOffsetStdev,
        noise_entries: t.noiseEntries.map((entry) => omit(entry, 'id')),
    };
    if (timingObj.noise_entries?.length === 0) {
        delete (timingObj as Partial<typeof timingObj>).noise_entries;
    }
    return cleanObject(timingObj);
};

export const isFileBackedAntennaPendingFile = (a: Antenna): boolean =>
    (a.pattern === 'xml' || a.pattern === 'file') &&
    (a.filename ?? '').trim().length === 0;

export const serializeAntenna = (a: Antenna): SerializedAntenna => {
    if (isFileBackedAntennaPendingFile(a)) {
        return cleanObject({
            id: a.id,
            name: a.name,
            pattern: 'isotropic',
            efficiency: a.efficiency,
        });
    }

    const rest = omit(a, 'type', 'meshScale');
    return cleanObject(rest) as SerializedAntenna;
};

export function deriveUtmZone(longitude: number): number {
    if (!Number.isFinite(longitude)) {
        return 1;
    }
    const clampedLongitude = Math.max(-180, Math.min(180, longitude));
    return Math.max(
        1,
        Math.min(60, Math.floor((clampedLongitude + 180) / 6) + 1)
    );
}

const serializeCoordinateSystem = (
    gp: GlobalParameters
): GlobalParameters['coordinateSystem'] => {
    if (gp.coordinateSystem.frame !== 'UTM') {
        return {
            frame: gp.coordinateSystem.frame,
        };
    }

    return {
        frame: 'UTM',
        zone: gp.coordinateSystem.zone ?? deriveUtmZone(gp.origin.longitude),
        hemisphere:
            gp.coordinateSystem.hemisphere ??
            (gp.origin.latitude < 0 ? 'S' : 'N'),
    };
};

export const serializeGlobalParameters = (gp: GlobalParameters) => {
    const { start, end, random_seed, oversample_ratio } = gp;
    const gpRest = omit(
        gp,
        'start',
        'end',
        'random_seed',
        'oversample_ratio',
        'coordinateSystem'
    );

    return cleanObject({
        ...gpRest,
        starttime: start,
        endtime: end,
        randomseed: random_seed,
        oversample: oversample_ratio,
        rotationangleunit: gp.rotationAngleUnit,
        coordinatesystem: serializeCoordinateSystem(gp),
    });
};
