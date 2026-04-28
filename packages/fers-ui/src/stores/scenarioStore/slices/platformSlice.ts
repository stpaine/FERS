// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { invoke } from '@tauri-apps/api/core';
import { StateCreator } from 'zustand';
import { createDefaultPlatform } from '../defaults';
import { generateSimId } from '../idUtils';
import { createUniqueScenarioName } from '../nameUtils';
import {
    cleanObject,
    serializeComponentInner,
    serializePlatform,
} from '../serializers';
import { enqueueFullSync, enqueueGranularSyncDetached } from '../syncQueue';
import {
    Platform,
    PlatformActions,
    PlatformComponent,
    ScenarioStore,
} from '../types';
import { buildScenarioJson } from './backendSlice';

const NUM_PATH_POINTS = 100;
type InterpolationType = 'static' | 'linear' | 'cubic';
interface InterpolatedPoint {
    x: number;
    y: number;
    z: number;
    vx: number;
    vy: number;
    vz: number;
}

interface InterpolatedRotationPoint {
    azimuth: number;
    elevation: number;
}

function syncPlatformGranular(
    getFn: () => ScenarioStore,
    platformId: string
): void {
    const platform = getFn().platforms.find((p) => p.id === platformId);
    if (!platform) return;
    enqueueGranularSyncDetached(
        'Platform',
        platform.id,
        JSON.stringify(cleanObject(serializePlatform(platform)))
    );
}

export const createPlatformSlice: StateCreator<
    ScenarioStore,
    [['zustand/immer', never]],
    [],
    PlatformActions
> = (set, get) => ({
    addPlatform: () => {
        set((state) => {
            const id = generateSimId('Platform');
            const newName = `Platform ${state.platforms.length + 1}`;
            const newPlatform: Platform = {
                ...createDefaultPlatform(),
                id,
                name: createUniqueScenarioName(state, newName),
            };
            // Defaults to empty components list
            state.platforms.push(newPlatform);
            state.isDirty = true;
        });
        // libfers has no granular add API for Platforms — full sync is required.
        void enqueueFullSync(() => buildScenarioJson(get()));
    },
    addPositionWaypoint: (platformId) => {
        let touched = false;
        set((state) => {
            const platform = state.platforms.find((p) => p.id === platformId);
            if (platform) {
                platform.motionPath.waypoints.push({
                    id: generateSimId('Platform'),
                    x: 0,
                    y: 0,
                    altitude: 0,
                    time: 0,
                });
                state.isDirty = true;
                touched = true;
            }
        });
        if (touched) syncPlatformGranular(get, platformId);
    },
    removePositionWaypoint: (platformId, waypointId) => {
        let touched = false;
        set((state) => {
            const platform = state.platforms.find((p) => p.id === platformId);
            if (platform && platform.motionPath.waypoints.length > 1) {
                const index = platform.motionPath.waypoints.findIndex(
                    (wp) => wp.id === waypointId
                );
                if (index > -1) {
                    platform.motionPath.waypoints.splice(index, 1);
                    state.isDirty = true;
                    touched = true;
                }
            }
        });
        if (touched) syncPlatformGranular(get, platformId);
    },
    addRotationWaypoint: (platformId) => {
        let touched = false;
        set((state) => {
            const platform = state.platforms.find((p) => p.id === platformId);
            if (platform?.rotation.type === 'path') {
                platform.rotation.waypoints.push({
                    id: generateSimId('Platform'),
                    azimuth: 0,
                    elevation: 0,
                    time: 0,
                });
                state.isDirty = true;
                touched = true;
            }
        });
        if (touched) syncPlatformGranular(get, platformId);
    },
    removeRotationWaypoint: (platformId, waypointId) => {
        let touched = false;
        set((state) => {
            const platform = state.platforms.find((p) => p.id === platformId);
            if (
                platform?.rotation.type === 'path' &&
                platform.rotation.waypoints.length > 1
            ) {
                const index = platform.rotation.waypoints.findIndex(
                    (wp) => wp.id === waypointId
                );
                if (index > -1) {
                    platform.rotation.waypoints.splice(index, 1);
                    state.isDirty = true;
                    touched = true;
                }
            }
        });
        if (touched) syncPlatformGranular(get, platformId);
    },
    addPlatformComponent: (platformId, componentType) => {
        let added = false;
        set((state) => {
            const platform = state.platforms.find((p) => p.id === platformId);
            if (!platform) return;

            const id =
                componentType === 'transmitter'
                    ? generateSimId('Transmitter')
                    : componentType === 'receiver'
                      ? generateSimId('Receiver')
                      : componentType === 'target'
                        ? generateSimId('Target')
                        : generateSimId('Transmitter');
            const rxId =
                componentType === 'monostatic'
                    ? generateSimId('Receiver')
                    : null;
            const name = `${platform.name} ${
                componentType.charAt(0).toUpperCase() + componentType.slice(1)
            }`;
            const uniqueName = createUniqueScenarioName(state, name);
            let newComponent: PlatformComponent;

            switch (componentType) {
                case 'monostatic':
                    newComponent = {
                        id,
                        type: 'monostatic',
                        txId: id,
                        rxId: rxId ?? generateSimId('Receiver'),
                        name: uniqueName,
                        radarType: 'pulsed',
                        window_skip: 0,
                        window_length: 1e-5,
                        prf: 1000,
                        antennaId: null,
                        waveformId: null,
                        timingId: null,
                        noiseTemperature: 290,
                        noDirectPaths: false,
                        noPropagationLoss: false,
                        schedule: [],
                    };
                    break;
                case 'transmitter':
                    newComponent = {
                        id,
                        type: 'transmitter',
                        name: uniqueName,
                        radarType: 'pulsed',
                        prf: 1000,
                        antennaId: null,
                        waveformId: null,
                        timingId: null,
                        schedule: [],
                    };
                    break;
                case 'receiver':
                    newComponent = {
                        id,
                        type: 'receiver',
                        name: uniqueName,
                        radarType: 'pulsed',
                        window_skip: 0,
                        window_length: 1e-5,
                        prf: 1000,
                        antennaId: null,
                        timingId: null,
                        noiseTemperature: 290,
                        noDirectPaths: false,
                        noPropagationLoss: false,
                        schedule: [],
                    };
                    break;
                case 'target':
                    newComponent = {
                        id,
                        type: 'target',
                        name: uniqueName,
                        rcs_type: 'isotropic',
                        rcs_value: 1,
                        rcs_model: 'constant',
                    };
                    break;
                default:
                    return;
            }
            platform.components.push(newComponent);
            state.isDirty = true;
            added = true;
        });
        if (added) {
            // Components (Transmitter/Receiver/Target/Monostatic) are independent
            // backend objects with their own IDs; there is no granular add API.
            void enqueueFullSync(() => buildScenarioJson(get()));
        }
    },
    removePlatformComponent: (platformId, componentId) => {
        let removed = false;
        set((state) => {
            const platform = state.platforms.find((p) => p.id === platformId);
            if (platform) {
                const index = platform.components.findIndex(
                    (c) => c.id === componentId
                );
                if (index > -1) {
                    platform.components.splice(index, 1);
                    if (state.selectedComponentId === componentId) {
                        state.selectedComponentId = null;
                    }
                    state.isDirty = true;
                    removed = true;
                }
            }
        });
        if (removed) {
            // libfers has no granular remove API — full sync is required.
            void enqueueFullSync(() => buildScenarioJson(get()));
        }
    },
    setPlatformRcsModel: (platformId, componentId, newModel) => {
        let touched = false;
        set((state) => {
            const platform = state.platforms.find((p) => p.id === platformId);
            const component = platform?.components.find(
                (c) => c.id === componentId
            );
            if (component?.type === 'target') {
                component.rcs_model = newModel;
                if (newModel === 'chisquare' || newModel === 'gamma') {
                    if (typeof component.rcs_k !== 'number') {
                        component.rcs_k = 1.0;
                    }
                } else {
                    delete component.rcs_k;
                }
                state.isDirty = true;
                touched = true;
            }
        });
        if (touched) {
            const platform = get().platforms.find((p) => p.id === platformId);
            const component = platform?.components.find(
                (c) => c.id === componentId
            );
            if (component) {
                enqueueGranularSyncDetached(
                    'Target',
                    component.id,
                    JSON.stringify(
                        cleanObject(serializeComponentInner(component))
                    )
                );
            }
        }
    },
    fetchPlatformPath: async (platformId) => {
        const { platforms, showError, globalParameters } = get();
        const platform = platforms.find((p) => p.id === platformId);

        if (!platform) return;

        // 1. Fetch/Calculate Motion Path
        const { waypoints, interpolation } = platform.motionPath;
        let newPathPoints: {
            x: number;
            y: number;
            z: number;
            vx: number;
            vy: number;
            vz: number;
        }[] = [];

        try {
            if (waypoints.length < 2 || interpolation === 'static') {
                // Static or single point: Calculate directly on frontend (velocity 0)
                newPathPoints = waypoints.map((wp) => ({
                    x: wp.x,
                    y: wp.altitude,
                    z: -wp.y, // ENU Y -> Three JS -Z
                    vx: 0,
                    vy: 0,
                    vz: 0,
                }));
            } else {
                // Dynamic: Fetch interpolated points from Backend
                const points = await invoke<InterpolatedPoint[]>(
                    'get_interpolated_motion_path',
                    {
                        waypoints,
                        interpType: interpolation as InterpolationType,
                        numPoints: NUM_PATH_POINTS,
                    }
                );
                // Convert ENU (Backend) to Three.js coordinates
                // Pos: X->X, Alt->Y, Y->-Z
                // Vel: Vx->Vx, Vz->Vy, Vy->-Vz
                newPathPoints = points.map((p) => ({
                    x: p.x,
                    y: p.z,
                    z: -p.y,
                    vx: p.vx,
                    vy: p.vz,
                    vz: -p.vy,
                }));
            }
        } catch (error) {
            const msg = error instanceof Error ? error.message : String(error);
            console.error(
                `Failed to fetch motion path for ${platform.name}:`,
                msg
            );
            showError(`Failed to get motion path for ${platform.name}: ${msg}`);
            // Fallback to empty to prevent stale data
            newPathPoints = [];
        }

        // 2. Fetch/Calculate Rotation Path
        const { rotation } = platform;
        let newRotationPoints:
            | { azimuth: number; elevation: number }[]
            | undefined = undefined;

        if (rotation.type === 'path') {
            const rotWaypoints = rotation.waypoints;
            // Only fetch from backend if dynamic. Static/Single points are handled by the real-time calculator.
            if (
                rotWaypoints.length >= 2 &&
                rotation.interpolation !== 'static'
            ) {
                try {
                    const points = await invoke<InterpolatedRotationPoint[]>(
                        'get_interpolated_rotation_path',
                        {
                            waypoints: rotWaypoints,
                            interpType:
                                rotation.interpolation as InterpolationType,
                            angleUnit: globalParameters.rotationAngleUnit,
                            numPoints: NUM_PATH_POINTS,
                        }
                    );
                    newRotationPoints = points.map((p) => ({
                        azimuth: p.azimuth,
                        elevation: p.elevation,
                    }));
                } catch (error) {
                    const msg =
                        error instanceof Error ? error.message : String(error);
                    console.error(
                        `Failed to fetch rotation path for ${platform.name}:`,
                        msg
                    );
                    // Log error but don't break the whole update; standard calc will fallback to first waypoint
                }
            }
        }

        // 3. Update Store
        set((state) => {
            const p = state.platforms.find((p) => p.id === platformId);
            if (p) {
                p.pathPoints = newPathPoints;
                p.rotationPathPoints = newRotationPoints;
            }
        });
    },
});
