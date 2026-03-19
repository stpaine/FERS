// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { useMemo, useRef } from 'react';
import * as THREE from 'three';
import { fersColors } from '@/theme';
import { useDynamicScale } from '@/hooks/useDynamicScale';

const BORESIGHT_LENGTH = 7; // Length of the arrow in world units.

/**
 * Renders a boresight arrow indicating the forward direction of an object.
 * This component should be placed inside a <group> that has its rotation set
 * to the object's orientation. The arrow points along the group's local -Z axis.
 */
export function BoresightArrow() {
    const groupRef = useRef<THREE.Group>(null);

    // Apply dynamic scaling
    useDynamicScale(groupRef);

    const arrowHelper = useMemo(() => {
        const dir = new THREE.Vector3(0, 0, -1); // Points "forward" along local -Z axis.
        const origin = new THREE.Vector3(0, 0, 0); // Arrow starts at the object's origin.
        return new THREE.ArrowHelper(
            dir,
            origin,
            BORESIGHT_LENGTH,
            fersColors.physics.boresight
        );
    }, []);

    // Use the <primitive> element to add the pre-built ArrowHelper object to the scene.
    return (
        <group ref={groupRef}>
            <primitive object={arrowHelper} />
        </group>
    );
}
