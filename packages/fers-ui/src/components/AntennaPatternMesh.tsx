// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { invoke } from '@tauri-apps/api/core';
import { useEffect, useMemo, useRef, useState } from 'react';
import * as THREE from 'three';
import { useShallow } from 'zustand/react/shallow';
import { useDynamicScale } from '@/hooks/useDynamicScale';
import { PlatformComponent, useScenarioStore } from '@/stores/scenarioStore';
import { isFileBackedAntennaPendingFile } from '@/stores/scenarioStore/serializers';
import { useSimulationProgressStore } from '@/stores/simulationProgressStore';

const AZIMUTH_SEGMENTS = 64; // Resolution for azimuth sampling
const ELEVATION_SEGMENTS = 32; // Resolution for elevation sampling
const BASE_MESH_RADIUS = 3; // Base visual scaling factor for the main lobe.
const PREVIEW_RETRY_DELAY_MS = 250;

export const BACKEND_BUSY_MESSAGE = 'Backend is busy with another operation';

interface AntennaPatternData {
    gains: number[];
    az_count: number;
    el_count: number;
    max_gain: number;
}

interface AntennaPatternMeshProps {
    antennaId: string;
    component: PlatformComponent;
}

type AntennaPreviewBusyState = {
    hasRequest: boolean;
    isBackendSyncing: boolean;
    isSimulating: boolean;
    isGeneratingKml: boolean;
};

type AntennaPreviewErrorAction = {
    clearError: boolean;
    clearPattern: boolean;
    scheduleRetry: boolean;
};

export function isBackendBusyError(error: unknown) {
    const message = error instanceof Error ? error.message : String(error);
    return message.includes(BACKEND_BUSY_MESSAGE);
}

export function getAntennaPreviewErrorAction(
    error: unknown
): AntennaPreviewErrorAction {
    const isBusy = isBackendBusyError(error);
    return {
        clearError: isBusy,
        clearPattern: !isBusy,
        scheduleRetry: isBusy,
    };
}

export function shouldDeferAntennaPreviewFetch({
    hasRequest,
    isBackendSyncing,
    isSimulating,
    isGeneratingKml,
}: AntennaPreviewBusyState) {
    return hasRequest && (isBackendSyncing || isSimulating || isGeneratingKml);
}

/**
 * Creates and renders a 3D heatmap mesh from real antenna pattern data
 * fetched from the `libfers` backend.
 *
 * @returns A Three.js mesh component representing the antenna gain pattern.
 */
export function AntennaPatternMesh({
    antennaId,
    component,
}: AntennaPatternMeshProps) {
    const [patternData, setPatternData] = useState<AntennaPatternData | null>(
        null
    );

    const isSimulating = useSimulationProgressStore(
        (state) => state.isSimulating
    );
    const isGeneratingKml = useSimulationProgressStore(
        (state) => state.isGeneratingKml
    );

    // Select antenna, potential waveform, and backend sync state
    const { antenna, waveform, isBackendSyncing } = useScenarioStore(
        useShallow((state) => {
            const ant = state.antennas.find((a) => a.id === antennaId);
            const wf =
                'waveformId' in component && component.waveformId
                    ? state.waveforms.find((w) => w.id === component.waveformId)
                    : undefined;
            return {
                antenna: ant,
                waveform: wf,
                isBackendSyncing: state.isBackendSyncing,
            };
        })
    );
    const { setAntennaPreviewError, clearAntennaPreviewError } =
        useScenarioStore.getState();

    const antennaIdStr = antenna?.id;
    const userScale = antenna?.meshScale ?? 1.0;
    const groupRef = useRef<THREE.Group>(null!);
    const retryTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);
    const lastRequestKeyRef = useRef<string | null>(null);
    const [retryNonce, setRetryNonce] = useState(0);

    // Create a stable hash of the antenna's properties to trigger real-time updates
    const antennaHash = JSON.stringify(antenna);

    // Apply dynamic scaling hook
    useDynamicScale(groupRef, { baseScale: userScale });

    // Determine the frequency to use for pattern calculation
    const frequency = useMemo(() => {
        if (!antenna) return null;
        if (isFileBackedAntennaPendingFile(antenna)) return null;

        // 1. If the antenna is frequency-independent, we can render it regardless of waveform.
        const independentTypes = [
            'isotropic',
            'sinc',
            'gaussian',
            'xml',
            'file',
        ];
        if (independentTypes.includes(antenna.pattern)) {
            // Use a dummy frequency (1GHz) if the backend requires non-zero,
            // effectively ignored by these pattern types.
            return 1e9;
        }

        // 2. If the antenna is frequency-dependent (Horn/Parabolic):
        // Priority A: Use the active Waveform frequency if attached.
        if (waveform?.carrier_frequency) {
            return waveform.carrier_frequency;
        }

        // Priority B: Use the 'Design Frequency' from the Antenna Inspector.
        if (antenna.design_frequency) {
            return antenna.design_frequency;
        }

        // 3. If prerequisites are not met, return null to suppress rendering.
        return null;
    }, [antenna, waveform]);
    const requestKey =
        antennaIdStr && frequency !== null
            ? `${antennaIdStr}:${frequency}:${antennaHash}`
            : null;
    const shouldDeferFetch = shouldDeferAntennaPreviewFetch({
        hasRequest: requestKey !== null,
        isBackendSyncing,
        isSimulating,
        isGeneratingKml,
    });

    useEffect(() => {
        return () => {
            if (retryTimerRef.current !== null) {
                clearTimeout(retryTimerRef.current);
                retryTimerRef.current = null;
            }
        };
    }, []);

    useEffect(() => {
        let isCancelled = false;

        if (retryTimerRef.current !== null) {
            clearTimeout(retryTimerRef.current);
            retryTimerRef.current = null;
        }

        if (!requestKey || !antennaIdStr || frequency === null) {
            lastRequestKeyRef.current = null;
            setPatternData(null);
            if (antennaIdStr) {
                clearAntennaPreviewError(antennaIdStr);
            }
            return () => {
                isCancelled = true;
            };
        }

        const requestKeyChanged = lastRequestKeyRef.current !== requestKey;
        lastRequestKeyRef.current = requestKey;
        if (requestKeyChanged) {
            setPatternData(null);
        }

        if (shouldDeferFetch) {
            clearAntennaPreviewError(antennaIdStr);
            return () => {
                isCancelled = true;
            };
        }

        const fetchPattern = async () => {
            try {
                const data = await invoke<AntennaPatternData>(
                    'get_antenna_pattern',
                    {
                        antennaId: antennaIdStr,
                        azSamples: AZIMUTH_SEGMENTS + 1,
                        elSamples: ELEVATION_SEGMENTS + 1,
                        frequency,
                    }
                );

                if (!isCancelled) {
                    clearAntennaPreviewError(antennaIdStr);
                    setPatternData(data);
                }
            } catch (error) {
                if (isCancelled) {
                    return;
                }

                const message =
                    error instanceof Error ? error.message : String(error);
                const action = getAntennaPreviewErrorAction(error);

                if (action.clearError) {
                    clearAntennaPreviewError(antennaIdStr);
                }

                if (action.clearPattern) {
                    console.error(
                        `Failed to fetch pattern for antenna ${antennaIdStr}:`,
                        error
                    );
                    setAntennaPreviewError(antennaIdStr, message);
                    setPatternData(null);
                }

                if (action.scheduleRetry) {
                    retryTimerRef.current = setTimeout(() => {
                        setRetryNonce((current) => current + 1);
                    }, PREVIEW_RETRY_DELAY_MS);
                }
            }
        };

        void fetchPattern();

        return () => {
            isCancelled = true;
        };
    }, [antennaIdStr, frequency, requestKey, retryNonce, shouldDeferFetch]);

    const geometry = useMemo(() => {
        if (!patternData) return new THREE.BufferGeometry();

        const geom = new THREE.BufferGeometry();
        const vertices: number[] = [];
        const colors: number[] = [];
        const { gains, az_count, el_count } = patternData;

        for (let i = 0; i < el_count; i++) {
            // Elevation from -PI/2 to PI/2
            const elevation = (i / (el_count - 1)) * Math.PI - Math.PI / 2;

            for (let j = 0; j < az_count; j++) {
                // Azimuth from -PI to PI
                const azimuth = (j / (az_count - 1)) * 2 * Math.PI - Math.PI;

                const gain = gains[i * az_count + j]; // Normalized gain [0, 1]
                const radius = gain * BASE_MESH_RADIUS;

                // Convert spherical to Cartesian for a -Z forward orientation
                const x = radius * Math.cos(elevation) * Math.sin(azimuth);
                const y = radius * Math.sin(elevation);
                const z = -radius * Math.cos(elevation) * Math.cos(azimuth);

                vertices.push(x, y, z);

                // Assign color based on gain (heatmap: blue -> green -> red)
                const color = new THREE.Color();
                // HSL: Hue from blue (0.66, low gain) to red (0.0, high gain).
                color.setHSL(0.66 * (1 - gain), 1.0, 0.5);
                colors.push(color.r, color.g, color.b);
            }
        }

        const indices: number[] = [];
        for (let i = 0; i < ELEVATION_SEGMENTS; i++) {
            for (let j = 0; j < AZIMUTH_SEGMENTS; j++) {
                const a = i * (AZIMUTH_SEGMENTS + 1) + j;
                const b = a + AZIMUTH_SEGMENTS + 1;
                const c = a + 1;
                const d = b + 1;

                indices.push(a, b, c); // First triangle
                indices.push(b, d, c); // Second triangle
            }
        }

        geom.setIndex(indices);
        geom.setAttribute(
            'position',
            new THREE.Float32BufferAttribute(vertices, 3)
        );
        geom.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        geom.computeVertexNormals();

        return geom;
    }, [patternData]);

    if (!patternData) return null;

    return (
        <group ref={groupRef}>
            <mesh geometry={geometry}>
                <meshStandardMaterial
                    vertexColors
                    transparent
                    opacity={0.5}
                    side={THREE.DoubleSide}
                    roughness={0.7}
                    metalness={0.1}
                    depthWrite={false}
                    polygonOffset
                    polygonOffsetFactor={1}
                    polygonOffsetUnits={1}
                />
            </mesh>
        </group>
    );
}
