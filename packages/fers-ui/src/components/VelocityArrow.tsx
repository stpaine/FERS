// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { useMemo, useRef } from 'react';
import * as THREE from 'three';
import {
    Platform,
    calculateInterpolatedVelocity,
} from '@/stores/scenarioStore';
import { fersColors } from '@/theme';
import { useDynamicScale } from '@/hooks/useDynamicScale';

const VELOCITY_ARROW_LENGTH = 7; // Fixed base length (matches Boresight)

interface VelocityArrowProps {
    platform: Platform;
    currentTime: number;
}

/**
 * Renders an arrow indicating the velocity vector of the platform.
 * The direction indicates movement direction.
 * Note: Length is fixed for visualization; actual speed is in properties.
 */
export function VelocityArrow({ platform, currentTime }: VelocityArrowProps) {
    const groupRef = useRef<THREE.Group>(null);

    // Apply dynamic scaling
    useDynamicScale(groupRef);

    const velocity = useMemo(
        () => calculateInterpolatedVelocity(platform, currentTime),
        [platform, currentTime]
    );

    const arrowHelper = useMemo(() => {
        const speed = velocity.length();
        // Don't render if static
        if (speed < 0.001) return null;

        const dir = velocity.clone().normalize();
        const origin = new THREE.Vector3(0, 0, 0); // Local origin of the container

        return new THREE.ArrowHelper(
            dir,
            origin,
            VELOCITY_ARROW_LENGTH,
            fersColors.physics.velocity
        );
    }, [velocity]);

    if (!arrowHelper) return null;

    return (
        <group ref={groupRef}>
            <primitive object={arrowHelper} />
        </group>
    );
}
