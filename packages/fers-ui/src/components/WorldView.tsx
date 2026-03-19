// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import React, { useEffect, useMemo, useRef, memo } from 'react';
import { MapControls, Environment, Html } from '@react-three/drei';
import { Vector3 } from 'three';
import * as THREE from 'three';
import {
    useScenarioStore,
    Platform,
    calculateInterpolatedPosition,
    calculateInterpolatedRotation,
} from '@/stores/scenarioStore';
import { MotionPathLine } from './MotionPathLine';
import CameraManager from './CameraManager';
import { type MapControls as MapControlsImpl } from 'three-stdlib';
import { BoresightArrow } from './BoresightArrow';
import { VelocityArrow } from './VelocityArrow';
import { AntennaPatternMesh } from './AntennaPatternMesh';
import LinkVisualizer from './LinkVisualizer';
import { fersColors } from '@/theme';
import { useDynamicScale } from '@/hooks/useDynamicScale';

/**
 * A wrapper for axesHelper that scales based on camera distance to maintain visibility.
 */
function ScaledAxesHelper({
    visible,
    size = 2,
}: {
    visible: boolean;
    size?: number;
}) {
    const ref = useRef<THREE.Group>(null);

    // Apply dynamic scaling hook
    useDynamicScale(ref);

    return (
        <group ref={ref}>
            <axesHelper args={[size]} visible={visible} />
        </group>
    );
}

/**
 * Custom hook to calculate a platform's position at a given simulation time.
 * It relies on the pre-fetched path data stored in the Zustand store.
 * @param {Platform} platform The platform data.
 * @param {number} currentTime The current global simulation time.
 * @returns {Vector3} The interpolated position in Three.js coordinates.
 */
function useInterpolatedPosition(
    platform: Platform,
    currentTime: number
): Vector3 {
    return useMemo(
        () => calculateInterpolatedPosition(platform, currentTime),
        [platform, currentTime]
    );
}

/**
 * Custom hook to calculate a platform's rotation at a given simulation time.
 */
function useInterpolatedRotation(platform: Platform, currentTime: number) {
    return useMemo(
        () => calculateInterpolatedRotation(platform, currentTime),
        [platform, currentTime]
    );
}

/**
 * Represents a single platform in the 3D world as a sphere with an identifying label.
 * @param {object} props - The component props.
 * @param {Platform} props.platform - The platform data from the store.
 */
const PlatformSphere = memo(function PlatformSphere({
    platform,
}: {
    platform: Platform;
}) {
    const currentTime = useScenarioStore((state) => state.currentTime);

    const isSelected = useScenarioStore(
        (state) => state.selectedItemId === platform.id
    );

    // Visibility Toggles
    const {
        showPlatforms,
        showPlatformLabels,
        showAxes,
        showBoresights,
        showPatterns,
        showVelocities,
    } = useScenarioStore((state) => state.visibility);

    // Hooks must run unconditionally
    const position = useInterpolatedPosition(platform, currentTime);
    const rotation = useInterpolatedRotation(platform, currentTime);

    // Find all components on this platform that have an antenna
    const antennaComponents = useMemo(() => {
        return platform.components.filter(
            (
                c
            ): c is Extract<
                typeof c,
                { type: 'monostatic' | 'transmitter' | 'receiver' }
            > =>
                (c.type === 'monostatic' ||
                    c.type === 'transmitter' ||
                    c.type === 'receiver') &&
                c.antennaId !== null
        );
    }, [platform.components]);

    // Find all target components to visualize their RCS
    const targetComponents = useMemo(() => {
        return platform.components.filter(
            (c): c is Extract<typeof c, { type: 'target' }> =>
                c.type === 'target'
        );
    }, [platform.components]);

    const labelData = useMemo(
        () => ({
            x: position?.x,
            y: -position?.z, // Convert from Three.js Z back to ENU Y
            altitude: position?.y, // Convert from Three.js Y back to ENU Altitude
        }),
        [position]
    );

    return (
        <group position={position}>
            <group rotation={rotation}>
                <mesh visible={showPlatforms}>
                    <sphereGeometry args={[0.5, 32, 32]} />
                    <meshStandardMaterial
                        color={
                            isSelected
                                ? fersColors.platform.selected
                                : fersColors.platform.default
                        }
                        roughness={0.2}
                        metalness={0.8}
                        emissive={
                            isSelected
                                ? fersColors.platform.emission
                                : '#000000'
                        }
                        emissiveIntensity={isSelected ? 0.4 : 0}
                    />
                    {/* Render Body Axes: Red=X (Right), Green=Y (Up), Blue=Z (Rear). */}
                    <ScaledAxesHelper visible={showAxes} size={2} />
                </mesh>

                {/* Visualize Isotropic Static Sphere RCS */}
                {targetComponents.map((target) => {
                    {
                        /* TODO: currently only rendering constant isotropic RCS */
                    }
                    if (
                        target.rcs_type === 'isotropic' &&
                        target.rcs_value &&
                        target.rcs_value > 0
                    ) {
                        const radius = Math.sqrt(target.rcs_value / Math.PI);
                        return (
                            <mesh key={target.id} visible={showPlatforms}>
                                <sphereGeometry args={[radius, 24, 24]} />
                                <meshBasicMaterial
                                    color={fersColors.physics.rcs}
                                    wireframe
                                    transparent
                                    opacity={0.15}
                                />
                            </mesh>
                        );
                    }
                    return null;
                })}

                {/* Boresight is lightweight, conditional is fine, but visible is smoother */}
                <group
                    visible={
                        showPlatforms &&
                        showBoresights &&
                        antennaComponents.length > 0
                    }
                >
                    <BoresightArrow />
                </group>

                {/* Show antenna pattern meshes if the global toggle is enabled */}
                {showPlatforms &&
                    showPatterns &&
                    antennaComponents.map((comp) =>
                        comp.antennaId ? (
                            <AntennaPatternMesh
                                key={comp.id}
                                antennaId={comp.antennaId}
                                component={comp}
                            />
                        ) : null
                    )}
            </group>

            {/* Velocity Arrow */}
            <group visible={showPlatforms && showVelocities}>
                <VelocityArrow platform={platform} currentTime={currentTime} />
            </group>

            {/* Label */}
            {showPlatforms && showPlatformLabels && (
                <Html
                    position={[0, 1.2, 0]} // Position label above the sphere
                    center // Center the label on its anchor point
                    zIndexRange={[100, 0]} // Ensure labels render below UI overlays (zIndex 1000)
                    style={{
                        backgroundColor: fersColors.background.overlay,
                        color: fersColors.text.primary,
                        padding: '4px 8px',
                        borderRadius: '4px',
                        fontSize: '12px',
                        fontFamily: 'monospace',
                        whiteSpace: 'nowrap',
                        pointerEvents: 'none', // Allows camera clicks to pass through
                        userSelect: 'none', // Prevents text selection
                        transition: 'opacity 0.2s, border 0.2s',
                        opacity: isSelected ? 1 : 0.6,
                        border: `1px solid ${
                            isSelected
                                ? fersColors.platform.selected
                                : fersColors.text.secondary
                        }`,
                        boxShadow: isSelected
                            ? `0 0 8px ${fersColors.platform.selected}40`
                            : 'none',
                    }}
                >
                    <div style={{ fontWeight: 'bold' }}>{platform.name}</div>
                    <div>{`X: ${(labelData?.x ?? 0).toFixed(2)}`}</div>
                    <div>{`Y: ${(labelData?.y ?? 0).toFixed(2)}`}</div>
                    <div>{`Alt: ${(labelData?.altitude ?? 0).toFixed(2)}`}</div>
                </Html>
            )}
        </group>
    );
});

/**
 * WorldView represents the 3D scene where simulation elements are visualized.
 * It provides an interactive camera and renders platforms from the current scenario.
 */
interface WorldViewProps {
    controlsRef: React.RefObject<MapControlsImpl | null>;
}

/**
 * WorldView represents the 3D scene where simulation elements are visualized.
 * It provides an interactive camera and renders platforms from the current scenario.
 */
export default function WorldView({ controlsRef }: WorldViewProps) {
    const platforms = useScenarioStore((state) => state.platforms);

    const fetchPlatformPath = useScenarioStore(
        (state) => state.fetchPlatformPath
    );

    // Root Level Visibility Toggles
    const { showLinks, showMotionPaths } = useScenarioStore(
        (state) => state.visibility
    );

    // Keep a reference to platforms to access in the effect without adding it to dependencies.
    // We update this ref in a separate effect to avoid "mutation during render" errors.
    const platformsRef = useRef(platforms);

    useEffect(() => {
        platformsRef.current = platforms;
    }, [platforms]);

    // Memoize platform dependencies.
    const platformDeps = useMemo(
        () =>
            platforms
                .map((p) => {
                    const rotKey =
                        p.rotation.type === 'path'
                            ? `${p.rotation.interpolation}-${JSON.stringify(
                                  p.rotation.waypoints
                              )}`
                            : `fixed-${p.rotation.startAzimuth}-${p.rotation.startElevation}-${p.rotation.azimuthRate}-${p.rotation.elevationRate}`;

                    return [
                        p.id,
                        p.motionPath.interpolation,
                        JSON.stringify(p.motionPath.waypoints),
                        rotKey,
                    ].join('|');
                })
                .join(';'),
        [platforms]
    );

    useEffect(() => {
        // Access the platforms via ref. The dependency array is strictly controlled
        // by platformDeps, which prevents the infinite loop caused by store updates.
        platformsRef.current.forEach((platform) => {
            void fetchPlatformPath(platform.id);
        });
    }, [platformDeps, fetchPlatformPath]);

    return (
        <>
            <CameraManager controlsRef={controlsRef} />

            {/* Controls */}
            <MapControls makeDefault ref={controlsRef} />

            {/* Lighting */}
            <ambientLight intensity={0.5} />
            <directionalLight position={[50, 50, 25]} intensity={2.5} />

            {/* Environment for realistic reflections and ambient light */}
            <Environment files="/potsdamer_platz_1k.hdr" />

            {/* Physics Link Visualization - Conditional Render */}
            {showLinks && <LinkVisualizer />}

            {/* Objects */}
            {platforms.map((platform) => (
                <group key={platform.id}>
                    <PlatformSphere platform={platform} />
                    {/* Motion Path Lines */}
                    {showMotionPaths && <MotionPathLine platform={platform} />}
                </group>
            ))}
        </>
    );
}
