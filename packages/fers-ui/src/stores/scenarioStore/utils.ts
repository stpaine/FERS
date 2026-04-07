// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { Euler, Vector3 } from 'three';
import {
    GlobalParameters,
    Platform,
    PlatformComponent,
    ScenarioItem,
    ScenarioState,
} from './types';

// Helper to set nested properties safely
export const setPropertyByPath = (
    obj: object,
    path: string,
    value: unknown
): void => {
    const keys = path.split('.');
    const lastKey = keys.pop();
    if (!lastKey) return;

    let current: Record<string, unknown> = obj as Record<string, unknown>;

    for (const key of keys) {
        const next = current[key];
        if (typeof next !== 'object' || next === null) {
            // Path does not exist, so we cannot set the value.
            return;
        }
        current = next as Record<string, unknown>;
    }

    current[lastKey] = value;
};

// Helper function to find any item in the store by its ID
export const findItemInStore = (
    state: ScenarioState,
    id: string | null
): ScenarioItem | null => {
    if (!id) return null;
    if (id === 'global-parameters') return state.globalParameters;

    const collections = [
        state.waveforms,
        state.timings,
        state.antennas,
        state.platforms,
    ];
    for (const collection of collections) {
        const item = collection.find((i) => i.id === id);
        if (item) return item as ScenarioItem;
    }
    return null;
};

// Helper to find a platform component and its parent platform by component ID
export const findComponentInStore = (
    state: ScenarioState,
    componentId: string | null
): { platform: Platform; component: PlatformComponent } | null => {
    if (!componentId) return null;
    for (const platform of state.platforms) {
        const component = platform.components.find((c) => c.id === componentId);
        if (component) {
            return { platform, component };
        }
    }
    return null;
};

export function angleUnitToRadians(
    value: number,
    unit: GlobalParameters['rotationAngleUnit']
): number {
    return unit === 'deg' ? (value * Math.PI) / 180 : value;
}

/**
 * Calculates a platform's interpolated 3D position at a specific time.
 * This function relies on the pre-fetched `pathPoints` array stored on the platform object.
 * @param {Platform} platform The platform data, including its waypoints and cached path points.
 * @param {number} currentTime The global simulation time.
 * @returns {Vector3} The interpolated position in Three.js coordinates.
 */
export function calculateInterpolatedPosition(
    platform: Platform,
    currentTime: number
): Vector3 {
    const { waypoints, interpolation } = platform.motionPath;
    const pathPoints = platform.pathPoints ?? [];

    const firstWaypoint = waypoints[0];
    if (!firstWaypoint) return new Vector3(0, 0, 0);

    const staticPosition = new Vector3(
        firstWaypoint.x ?? 0,
        firstWaypoint.altitude ?? 0,
        -(firstWaypoint.y ?? 0)
    );

    if (
        interpolation === 'static' ||
        waypoints.length < 2 ||
        pathPoints.length < 2
    ) {
        return staticPosition;
    }

    const lastWaypoint = waypoints[waypoints.length - 1];
    const pathStartTime = firstWaypoint.time;
    const pathEndTime = lastWaypoint.time;
    const pathDuration = pathEndTime - pathStartTime;

    if (pathDuration <= 0) return staticPosition;

    const timeRatio = (currentTime - pathStartTime) / pathDuration;
    const clampedRatio = Math.max(0, Math.min(1, timeRatio));

    const floatIndex = clampedRatio * (pathPoints.length - 1);
    const index1 = Math.floor(floatIndex);
    const index2 = Math.min(pathPoints.length - 1, Math.ceil(floatIndex));

    const point1 = pathPoints[index1];
    const point2 = pathPoints[index2];

    if (!point1 || !point2) return staticPosition;

    const v1 = new Vector3(point1.x, point1.y, point1.z);
    if (index1 === index2) return v1;

    const v2 = new Vector3(point2.x, point2.y, point2.z);
    const interPointRatio = floatIndex - index1;
    return v1.clone().lerp(v2, interPointRatio);
}

/**
 * Calculates a platform's interpolated velocity vector at a specific time.
 * @param {Platform} platform The platform data.
 * @param {number} currentTime The global simulation time.
 * @returns {Vector3} The interpolated velocity in Three.js coordinates.
 */
export function calculateInterpolatedVelocity(
    platform: Platform,
    currentTime: number
): Vector3 {
    const { waypoints, interpolation } = platform.motionPath;
    const pathPoints = platform.pathPoints ?? [];
    const firstWaypoint = waypoints[0];

    if (
        !firstWaypoint ||
        interpolation === 'static' ||
        waypoints.length < 2 ||
        pathPoints.length < 2
    ) {
        return new Vector3(0, 0, 0);
    }

    const lastWaypoint = waypoints[waypoints.length - 1];
    const pathStartTime = firstWaypoint.time;
    const pathEndTime = lastWaypoint.time;
    const pathDuration = pathEndTime - pathStartTime;

    if (pathDuration <= 0) return new Vector3(0, 0, 0);

    const timeRatio = (currentTime - pathStartTime) / pathDuration;
    const clampedRatio = Math.max(0, Math.min(1, timeRatio));
    const floatIndex = clampedRatio * (pathPoints.length - 1);
    const index1 = Math.floor(floatIndex);
    const index2 = Math.min(pathPoints.length - 1, Math.ceil(floatIndex));

    const p1 = pathPoints[index1];
    const p2 = pathPoints[index2];

    if (!p1 || !p2) return new Vector3(0, 0, 0);
    if (index1 === index2) return new Vector3(p1.vx, p1.vy, p1.vz);

    const interPointRatio = floatIndex - index1;
    const v1 = new Vector3(p1.vx, p1.vy, p1.vz);
    const v2 = new Vector3(p2.vx, p2.vy, p2.vz);
    return v1.lerp(v2, interPointRatio);
}

/**
 * Calculates a platform's interpolated rotation (Euler) at a specific time.
 * @param {Platform} platform The platform data.
 * @param {number} currentTime The global simulation time.
 * @returns {Euler} The interpolated rotation in Three.js coordinates (YXZ order).
 */
export function calculateInterpolatedRotation(
    platform: Platform,
    currentTime: number,
    angleUnit: GlobalParameters['rotationAngleUnit']
): Euler {
    const { rotation } = platform;
    let azDeg = 0;
    let elDeg = 0;

    if (rotation.type === 'fixed') {
        // Linear calculation based on rate
        const dt = Math.max(0, currentTime); // Assume t=0 start for fixed
        azDeg = rotation.startAzimuth + rotation.azimuthRate * dt;
        elDeg = rotation.startElevation + rotation.elevationRate * dt;
    } else {
        // Path based interpolation
        const waypoints = rotation.waypoints;
        const pathPoints = platform.rotationPathPoints ?? [];
        const firstWaypoint = waypoints[0];

        if (!firstWaypoint) return new Euler(0, 0, 0);

        // Default to start
        azDeg = firstWaypoint.azimuth;
        elDeg = firstWaypoint.elevation;

        if (
            rotation.interpolation !== 'static' &&
            waypoints.length >= 2 &&
            pathPoints.length >= 2
        ) {
            const lastWaypoint = waypoints[waypoints.length - 1];
            const pathStartTime = firstWaypoint.time;
            const pathDuration = lastWaypoint.time - pathStartTime;

            if (pathDuration > 0) {
                const timeRatio = (currentTime - pathStartTime) / pathDuration;
                const clampedRatio = Math.max(0, Math.min(1, timeRatio));
                const floatIndex = clampedRatio * (pathPoints.length - 1);
                const index1 = Math.floor(floatIndex);
                const index2 = Math.min(
                    pathPoints.length - 1,
                    Math.ceil(floatIndex)
                );

                const p1 = pathPoints[index1];
                const p2 = pathPoints[index2];

                if (p1 && p2) {
                    const t = floatIndex - index1;
                    // Simple linear interpolation of angles for visualization
                    azDeg = p1.azimuth + (p2.azimuth - p1.azimuth) * t;
                    elDeg = p1.elevation + (p2.elevation - p1.elevation) * t;
                }
            }
        }
    }

    // Convert Compass Degrees (0 is North, CW) to Three.js Radians (0 is -Z?, CCW?)
    // FERS: 0 Az = North (Y), 90 Az = East (X).
    // Three.js: Y is Up.
    // We apply Azimuth as rotation around Y.
    // We apply Elevation as rotation around X.

    // Convert deg to rad
    const azRad = -angleUnitToRadians(azDeg, angleUnit); // Negate for CCW rotation in Three.js vs CW compass
    const elRad = angleUnitToRadians(elDeg, angleUnit);

    // Order YXZ: Rotate Azimuth (Y) first, then Elevation (X) (Pitch)
    return new Euler(elRad, azRad, 0, 'YXZ');
}
