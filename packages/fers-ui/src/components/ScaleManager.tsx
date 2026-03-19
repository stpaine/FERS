// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { useThree, useFrame } from '@react-three/fiber';
import { type MapControls as MapControlsImpl } from 'three-stdlib';
import * as THREE from 'three';
import React from 'react';

interface ScaleManagerProps {
    controlsRef: React.RefObject<MapControlsImpl | null>;
    labelRef: React.RefObject<HTMLDivElement | null>;
    barRef: React.RefObject<HTMLDivElement | null>;
}

/**
 * ScaleManager is a logic-only component that lives inside the Canvas.
 * It calculates the screen-to-world ratio and updates the external DOM elements
 * directly to ensure high performance without React renders.
 */
export function ScaleManager({
    controlsRef,
    labelRef,
    barRef,
}: ScaleManagerProps) {
    const { camera, gl } = useThree();

    // Target pixel width for the scale bar
    const targetWidthPx = 140;

    useFrame(() => {
        if (
            !controlsRef.current ||
            !labelRef.current ||
            !barRef.current ||
            !(camera instanceof THREE.PerspectiveCamera)
        ) {
            return;
        }

        const controls = controlsRef.current;

        // 1. Calculate distance from camera to the orbit target
        const distance = camera.position.distanceTo(controls.target);

        // 2. Calculate visible height at that distance (Vertical FOV)
        const vFOV = (camera.fov * Math.PI) / 180;
        const visibleHeightAtTarget = 2 * Math.tan(vFOV / 2) * distance;

        // 3. Calculate meters per pixel
        const canvasHeight = gl.domElement.clientHeight;
        const unitsPerPixel = visibleHeightAtTarget / canvasHeight;

        // 4. Calculate raw world units for our target pixel width
        const rawUnits = unitsPerPixel * targetWidthPx;

        // 5. Snap to "nice" numbers (1, 2, 5, 10...)
        const magnitude = Math.pow(10, Math.floor(Math.log10(rawUnits)));
        const residual = rawUnits / magnitude;

        let niceStep = 1;
        if (residual > 5) niceStep = 10;
        else if (residual > 2) niceStep = 5;
        else if (residual > 1) niceStep = 2;

        const niceUnits = niceStep * magnitude;

        // 6. Calculate exact pixel width for this nice unit amount
        const finalPixelWidth = niceUnits / unitsPerPixel;

        // 7. Format Label
        let labelText = '';
        if (niceUnits >= 1000) {
            labelText = `${(niceUnits / 1000).toFixed(0)} km`;
        } else {
            labelText = `${niceUnits.toFixed(0)} m`;
        }

        // 8. Direct DOM updates
        if (barRef.current.style.width !== `${finalPixelWidth}px`) {
            barRef.current.style.width = `${finalPixelWidth}px`;
        }
        if (labelRef.current.innerText !== labelText) {
            labelRef.current.innerText = labelText;
        }
    });

    return null;
}
