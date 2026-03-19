// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { useEffect, useRef } from 'react';
import { useThree, useFrame } from '@react-three/fiber';
import {
    useScenarioStore,
    calculateInterpolatedPosition,
} from '@/stores/scenarioStore';
import { type MapControls as MapControlsImpl } from 'three-stdlib';
import * as THREE from 'three';

interface CameraManagerProps {
    controlsRef: React.RefObject<MapControlsImpl | null>;
}

export default function CameraManager({ controlsRef }: CameraManagerProps) {
    const { camera } = useThree();

    const platforms = useScenarioStore((state) => state.platforms);
    const viewControlAction = useScenarioStore(
        (state) => state.viewControlAction
    );
    const clearViewControlAction = useScenarioStore(
        (state) => state.clearViewControlAction
    );

    const lastActionTimestamp = useRef(0);

    useEffect(() => {
        if (
            !controlsRef.current ||
            viewControlAction.timestamp === lastActionTimestamp.current
        ) {
            return;
        }

        const controls = controlsRef.current;
        const { type, targetId } = viewControlAction;

        if (type === 'frame') {
            const box = new THREE.Box3();
            let hasPoints = false;

            platforms.forEach((platform) => {
                (platform.pathPoints ?? []).forEach((point) => {
                    box.expandByPoint(
                        new THREE.Vector3(point.x, point.y, point.z)
                    );
                    hasPoints = true;
                });
            });

            if (
                hasPoints &&
                !box.isEmpty() &&
                camera instanceof THREE.PerspectiveCamera
            ) {
                const center = new THREE.Vector3();
                box.getCenter(center);
                const size = new THREE.Vector3();
                box.getSize(size);

                const maxDim = Math.max(size.x, size.y, size.z);
                const fov = camera.fov * (Math.PI / 180);
                let cameraZ = Math.abs(maxDim / 2 / Math.tan(fov / 2));
                cameraZ *= 1.5;

                const newPos = new THREE.Vector3(
                    center.x,
                    center.y + cameraZ / 2,
                    center.z + cameraZ
                );

                camera.position.copy(newPos);
                controls.target.copy(center);
                controls.update();
            }

            lastActionTimestamp.current = viewControlAction.timestamp;
            clearViewControlAction();
        } else if (type === 'focus' && targetId) {
            const platform = platforms.find((p) => p.id === targetId);
            if (platform) {
                const currentTime = useScenarioStore.getState().currentTime;
                const position = calculateInterpolatedPosition(
                    platform,
                    currentTime
                );
                controls.target.copy(position);
                controls.update();
            }
            lastActionTimestamp.current = viewControlAction.timestamp;
            clearViewControlAction();
        }
    }, [
        viewControlAction,
        controlsRef,
        camera,
        platforms,
        clearViewControlAction,
    ]);
    useFrame(() => {
        if (
            !controlsRef.current ||
            viewControlAction.type !== 'follow' ||
            !viewControlAction.targetId
        ) {
            return;
        }

        const controls = controlsRef.current;
        const platform = platforms.find(
            (p) => p.id === viewControlAction.targetId
        );

        if (platform) {
            const position = calculateInterpolatedPosition(
                platform,
                useScenarioStore.getState().currentTime
            );
            controls.target.copy(position);
            controls.update();
        }
    });

    return null;
}
