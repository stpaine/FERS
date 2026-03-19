// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { useEffect, useMemo, useState } from 'react';
import { Line } from '@react-three/drei';
import { Vector3 } from 'three';
import { invoke } from '@tauri-apps/api/core';
import { useScenarioStore, Platform } from '@/stores/scenarioStore';
import { fersColors } from '@/theme';

const NUM_PATH_POINTS = 100; // The resolution of the rendered path line.

// --- Type definitions for Tauri backend communication ---
type InterpolationType = 'static' | 'linear' | 'cubic';

interface InterpolatedPoint {
    x: number; // Corresponds to ENU X -> Three.js X
    y: number; // Corresponds to ENU Y -> Three.js -Z
    z: number; // Corresponds to ENU Altitude -> Three.js Y
}

/**
 * A component that fetches and renders the interpolated motion path for a given platform.
 * It communicates with the Tauri backend to compute the path using the libfers core.
 * @param {object} props - The component props.
 * @param {Platform} props.platform - The platform whose motion path should be rendered.
 */
export function MotionPathLine({ platform }: { platform: Platform }) {
    const [pathPoints, setPathPoints] = useState<Vector3[] | null>(null);
    const showError = useScenarioStore((state) => state.showError);

    const { waypoints, interpolation } = platform.motionPath;

    useEffect(() => {
        const fetchPath = async () => {
            if (waypoints.length === 0) {
                setPathPoints(null);
                return;
            }

            // For static or single-waypoint paths, no need to call the backend.
            // A static object has no visible path.
            if (interpolation === 'static' || waypoints.length < 2) {
                setPathPoints(null);
                return;
            }

            try {
                // The frontend 'PositionWaypoint' has an 'altitude' field, which the
                // Rust backend expects. The Rust `MotionWaypoint` struct will map this
                // correctly to the C-API's 'z' field.
                const points = await invoke<InterpolatedPoint[]>(
                    'get_interpolated_motion_path',
                    {
                        waypoints: waypoints,
                        interpType: interpolation as InterpolationType,
                        numPoints: NUM_PATH_POINTS,
                    }
                );

                const vectors = points.map(
                    // ENU to Three.js coordinate system mapping:
                    // Backend X (East) -> Three.js X
                    // Backend Z (Up) -> Three.js Y
                    // Backend Y (North) -> Three.js -Z
                    (p) => new Vector3(p.x, p.z, -p.y)
                );
                setPathPoints(vectors);
            } catch (error) {
                const errorMessage =
                    error instanceof Error ? error.message : String(error);
                console.error(
                    `Failed to fetch motion path for platform ${platform.name}:`,
                    errorMessage
                );
                showError(
                    `Failed to get motion path for ${platform.name}: ${errorMessage}`
                );
                setPathPoints(null);
            }
        };

        void fetchPath();
    }, [waypoints, interpolation, platform.name, showError]);

    const linePoints = useMemo(() => {
        if (!pathPoints || pathPoints.length < 2) return undefined;
        return pathPoints;
    }, [pathPoints]);

    if (!linePoints) {
        return null;
    }

    return (
        <Line
            points={linePoints}
            color={fersColors.physics.motionPath}
            lineWidth={1.5}
            dashed={false}
        />
    );
}
