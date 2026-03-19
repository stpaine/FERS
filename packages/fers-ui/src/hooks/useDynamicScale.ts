// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { useMemo } from 'react';
import { useFrame, useThree } from '@react-three/fiber';
import * as THREE from 'three';

interface DynamicScaleOptions {
    /**
     * The base multiplier for the scale.
     * Defaults to 1.0.
     */
    baseScale?: number;
    /**
     * The distance at which the object appears at its "natural" size (scale = baseScale).
     * Defaults to 50.
     */
    referenceDistance?: number;
}

/**
 * A hook that dynamically scales a Three.js object based on its distance from the camera.
 * This ensures the object maintains a relatively constant screen size (billboard-like scaling)
 * while preserving its 3D orientation.
 *
 * @param objectRef A React ref pointing to the Three.js Object3D (Group/Mesh).
 * @param options Configuration options for scaling behavior.
 */
export function useDynamicScale(
    objectRef: React.RefObject<THREE.Object3D | null>,
    options: DynamicScaleOptions = {}
) {
    const { baseScale = 1.0, referenceDistance = 50 } = options;
    const { camera } = useThree();

    // Reusable vector to prevent garbage collection churn in the render loop
    const worldPos = useMemo(() => new THREE.Vector3(), []);

    useFrame(() => {
        if (!objectRef.current) return;

        // Get absolute world position to calculate accurate distance to camera
        objectRef.current.getWorldPosition(worldPos);
        const distance = camera.position.distanceTo(worldPos);

        // Calculate scale factor: (Distance / RefDist) * Base
        // Clamp minimum to prevent mathematical singularities or invisibility
        const scale = Math.max(
            0.001,
            (distance / referenceDistance) * baseScale
        );

        objectRef.current.scale.set(scale, scale, scale);
    });
}
