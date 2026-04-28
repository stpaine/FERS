// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { Box, Tooltip, Typography } from '@mui/material';
import { Html, Line } from '@react-three/drei';
import { useFrame, useThree } from '@react-three/fiber';
import { invoke } from '@tauri-apps/api/core';
import React, { useEffect, useMemo, useRef, useState } from 'react';
import * as THREE from 'three';
import {
    calculateInterpolatedPosition,
    Platform,
    useScenarioStore,
} from '@/stores/scenarioStore';
import { useSimulationProgressStore } from '@/stores/simulationProgressStore';
import { fersColors } from '@/theme';

const TYPE_MAP = ['monostatic', 'illuminator', 'scattered', 'direct'] as const;
const QUALITY_MAP = ['strong', 'weak'] as const;
const COLORS = {
    monostatic: {
        // Strong monostatic return (received power above noise floor)
        strong: fersColors.link.monostatic.strong,
        // Weak monostatic return (received power below noise floor - rendered as transparent "ghost" line)
        weak: fersColors.link.monostatic.weak,
    },
    // Shows power density at the target in dBW/m²
    illuminator: fersColors.link.illuminator,
    // Shows received power in dBm (Rendered as transparent "ghost" line if signal is below noise floor)
    scattered: fersColors.link.scattered,
    // Shows direct interference/leakage power in dBm
    direct: fersColors.link.direct,
};

// Metadata only - no Vector3 positions here.
// Positions are derived at 60FPS during render.
interface LinkMetadata {
    link_type: 'monostatic' | 'illuminator' | 'scattered' | 'direct';
    quality: 'strong' | 'weak';
    label: string;
    source_id: string;
    dest_id: string;
    origin_id: string;
    rcs: number; // RCS in m^2; negative if not applicable
    actual_power_dbm: number; // Received power with actual RCS; -999 if not applicable
    display_value: number; // Numeric value represented by label, in the label's unit
}

interface LinkRange {
    powerMin: number;
    powerMax: number;
    rcsMin: number;
    rcsMax: number;
}

// Derived object used for actual rendering
interface RenderableLink extends LinkMetadata {
    start: THREE.Vector3;
    end: THREE.Vector3;
    distance: number;
    source_name: string;
    dest_name: string;
    origin_name: string;
    fmcwDutyCycle: number | null;
    range: LinkRange | null;
}

// Define the shape coming from Rust
interface RustVisualLink {
    link_type: number;
    quality: number;
    label: string;
    source_id: string;
    dest_id: string;
    origin_id: string;
    rcs: number; // RCS in m^2; negative if not applicable
    actual_power_dbm: number; // Received power with actual RCS; -999 if not applicable
    display_value: number; // Numeric value represented by label, in the label's unit
}

// Helper to determine color based on link type and quality
const getLinkColor = (link: LinkMetadata) => {
    if (link.link_type === 'monostatic') {
        return link.quality === 'strong'
            ? COLORS.monostatic.strong
            : COLORS.monostatic.weak;
    } else if (link.link_type === 'illuminator') {
        return COLORS.illuminator;
    } else if (link.link_type === 'scattered') {
        return COLORS.scattered;
    } else if (link.link_type === 'direct') {
        return COLORS.direct;
    }
    return '#ffffff';
};

const getLinkValueUnit = (link: LinkMetadata) =>
    link.link_type === 'illuminator' ? 'dBW/m²' : 'dBm';

function LinkLine({ link }: { link: RenderableLink }) {
    const points = useMemo(
        () => [link.start, link.end] as [THREE.Vector3, THREE.Vector3],
        [link.start, link.end]
    );

    const color = getLinkColor(link);

    // Ghost Line Logic:
    // If the link is geometrically valid but radiometrically weak (Sub-Noise),
    // render it transparently to indicate "possibility" without "detection".
    const isGhost =
        (link.link_type === 'monostatic' || link.link_type === 'scattered') &&
        link.quality === 'weak';

    const opacity = isGhost ? 0.1 : 1.0;

    return (
        <Line
            points={points}
            color={color}
            lineWidth={1.5}
            transparent
            opacity={opacity}
        />
    );
}

function LabelItem({ link, color }: { link: RenderableLink; color: string }) {
    const dutyCycle = link.fmcwDutyCycle ?? 1.0;
    const hasDutyCycledValue =
        link.fmcwDutyCycle !== null && dutyCycle > 0 && dutyCycle < 0.9995;
    const dutyCycleOffsetDb = hasDutyCycledValue
        ? -10.0 * Math.log10(dutyCycle)
        : 0.0;
    const peakDisplayValue =
        hasDutyCycledValue && link.display_value > -990
            ? link.display_value + dutyCycleOffsetDb
            : null;
    const actualPeakPowerDbm =
        hasDutyCycledValue && link.actual_power_dbm > -990
            ? link.actual_power_dbm + dutyCycleOffsetDb
            : null;

    return (
        <Tooltip
            arrow
            placement="top"
            title={
                <Box sx={{ p: 0.5 }}>
                    <Typography variant="subtitle2" sx={{ fontWeight: 'bold' }}>
                        Link Details
                    </Typography>
                    <Typography variant="caption" display="block">
                        <b>Path Segment:</b> {link.source_name} →{' '}
                        {link.dest_name}
                    </Typography>
                    {link.link_type === 'scattered' && (
                        <Typography variant="caption" display="block">
                            <b>Illuminator:</b> {link.origin_name}
                        </Typography>
                    )}
                    <Typography variant="caption" display="block">
                        <b>Type:</b> {link.link_type}
                    </Typography>
                    <Typography variant="caption" display="block">
                        <b>Distance:</b> {(link.distance / 1000).toFixed(2)} km
                    </Typography>
                    <Typography variant="caption" display="block">
                        <b>
                            {hasDutyCycledValue ? 'Average Value:' : 'Value:'}
                        </b>{' '}
                        {link.label}
                    </Typography>
                    {hasDutyCycledValue && (
                        <>
                            <Typography variant="caption" display="block">
                                <b>FMCW Duty Cycle:</b>{' '}
                                {(dutyCycle * 100).toFixed(1)}%
                            </Typography>
                            {peakDisplayValue !== null && (
                                <Typography
                                    variant="caption"
                                    display="block"
                                    sx={{ pl: 1.5, opacity: 0.7 }}
                                >
                                    peak equivalent:{' '}
                                    {peakDisplayValue.toFixed(1)}{' '}
                                    {getLinkValueUnit(link)}
                                </Typography>
                            )}
                        </>
                    )}
                    {link.actual_power_dbm > -990 && (
                        <>
                            <Typography variant="caption" display="block">
                                <b>
                                    {hasDutyCycledValue
                                        ? 'Actual Average Power:'
                                        : 'Actual Power:'}
                                </b>{' '}
                                {link.actual_power_dbm.toFixed(1)} dBm
                            </Typography>
                            {actualPeakPowerDbm !== null && (
                                <Typography
                                    variant="caption"
                                    display="block"
                                    sx={{ pl: 1.5, opacity: 0.7 }}
                                >
                                    actual peak: {actualPeakPowerDbm.toFixed(1)}{' '}
                                    dBm
                                </Typography>
                            )}
                            {link.range &&
                                link.range.powerMin < link.range.powerMax && (
                                    <Typography
                                        variant="caption"
                                        display="block"
                                        sx={{ pl: 1.5, opacity: 0.7 }}
                                    >
                                        range: {link.range.powerMin.toFixed(1)}{' '}
                                        to {link.range.powerMax.toFixed(1)} dBm
                                    </Typography>
                                )}
                        </>
                    )}
                    {link.rcs >= 0 && (
                        <>
                            <Typography variant="caption" display="block">
                                <b>RCS:</b> {link.rcs.toFixed(2)} m^2
                            </Typography>
                            {link.range &&
                                link.range.rcsMin < link.range.rcsMax && (
                                    <Typography
                                        variant="caption"
                                        display="block"
                                        sx={{ pl: 1.5, opacity: 0.7 }}
                                    >
                                        range: {link.range.rcsMin.toFixed(2)} to{' '}
                                        {link.range.rcsMax.toFixed(2)} m^2
                                    </Typography>
                                )}
                        </>
                    )}
                </Box>
            }
        >
            <div
                style={{
                    background: fersColors.background.overlay,
                    color: color,
                    padding: '2px 6px',
                    borderRadius: '4px',
                    fontSize: '10px',
                    fontFamily: 'monospace',
                    whiteSpace: 'nowrap',
                    border: `1px solid ${color}`,
                    boxShadow: `0 0 4px ${color}20`,
                    userSelect: 'none',
                    cursor: 'help',
                    pointerEvents: 'auto',
                }}
            >
                {link.label}
            </div>
        </Tooltip>
    );
}

// Component responsible for rendering a stack of labels at a specific 3D position
function LabelCluster({
    position,
    links,
    divRef,
}: {
    position: THREE.Vector3;
    links: RenderableLink[];
    divRef: React.Ref<HTMLDivElement>;
}) {
    return (
        <Html position={position} center zIndexRange={[100, 0]}>
            <div
                ref={divRef}
                style={{
                    display: 'flex',
                    flexDirection: 'column',
                    gap: '2px',
                    alignItems: 'center',
                    pointerEvents: 'none',
                }}
            >
                {links.map((link, i) => {
                    const color = getLinkColor(link);
                    return <LabelItem key={i} link={link} color={color} />;
                })}
            </div>
        </Html>
    );
}

// Minimum screen-space separation (px) enforced between label cluster centers.
const LABEL_OVERLAP_PX = 60;

/**
 * Renders radio frequency links between platforms.
 *
 * Performance Design:
 * 1. Metadata Fetching: 60 FPS. Labels and quality update at full rate.
 * 2. Range Tracking: Per-link min/max accumulated in a ref (no re-render cost); read in useMemo.
 * 3. Geometry Rendering: 60 FPS (driven by store.currentTime). Lines always stick to platforms.
 */
export default function LinkVisualizer() {
    const currentTime = useScenarioStore((state) => state.currentTime);
    const platforms = useScenarioStore((state) => state.platforms);
    const waveforms = useScenarioStore((state) => state.waveforms);
    const visibility = useScenarioStore((state) => state.visibility);
    const isSimulating = useSimulationProgressStore(
        (state) => state.isSimulating
    );
    const isGeneratingKml = useSimulationProgressStore(
        (state) => state.isGeneratingKml
    );
    const {
        showLinkLabels,
        showLinkMonostatic,
        showLinkIlluminator,
        showLinkScattered,
        showLinkDirect,
    } = visibility;

    const { camera, gl } = useThree();

    // Store only the metadata (connectivity/labels), not positions
    const [linkMetadata, setLinkMetadata] = useState<LinkMetadata[]>([]);

    // Accumulates observed min/max per link key — never causes re-renders, read in useMemo
    const linkRangesRef = useRef<Map<string, LinkRange>>(new Map());

    // Throttle control
    const lastFetchTimeRef = useRef<number>(0);
    const isFetchingRef = useRef<boolean>(false);
    const isMountedRef = useRef<boolean>(true);
    const FETCH_INTERVAL_MS = 16;

    // Track component mount status to prevent setting state on unmounted component
    useEffect(() => {
        isMountedRef.current = true;
        return () => {
            isMountedRef.current = false;
        };
    }, []);

    // Clear stale links when a backend operation holds the mutex
    useEffect(() => {
        if (isSimulating || isGeneratingKml) {
            setLinkMetadata([]);
            linkRangesRef.current.clear();
        }
    }, [isSimulating, isGeneratingKml]);

    // Build a lookup map: Component Name -> Parent Platform
    const componentToPlatform = useMemo(() => {
        const map = new Map<string, Platform>();
        platforms.forEach((p) => {
            // Map platform components (Tx, Rx, Tgt) to the platform
            p.components.forEach((c) => {
                map.set(c.id, p);
                if (c.type === 'monostatic') {
                    map.set(c.txId, p);
                    map.set(c.rxId, p);
                }
            });
        });
        return map;
    }, [platforms]);

    // Maps component/sub-component IDs to display names.
    // Monostatic txId/rxId are internal sub-IDs — map them to the parent component name.
    const componentToName = useMemo(() => {
        const map = new Map<string, string>();
        platforms.forEach((p) => {
            p.components.forEach((c) => {
                map.set(c.id, c.name);
                if (c.type === 'monostatic') {
                    map.set(c.txId, c.name);
                    map.set(c.rxId, c.name);
                }
            });
        });
        return map;
    }, [platforms]);

    const fmcwDutyCycleByTransmitter = useMemo(() => {
        const waveformById = new Map(waveforms.map((w) => [w.id, w]));
        const dutyCycleByTx = new Map<string, number>();

        platforms.forEach((p) => {
            p.components.forEach((c) => {
                if (c.type !== 'transmitter' && c.type !== 'monostatic') {
                    return;
                }

                const waveform = c.waveformId
                    ? waveformById.get(c.waveformId)
                    : undefined;
                if (!waveform) {
                    return;
                }

                const dutyCycle =
                    waveform.waveformType === 'fmcw_triangle'
                        ? 1
                        : waveform.waveformType === 'fmcw_linear_chirp' &&
                            waveform.chirp_period > 0
                          ? Math.min(
                                1,
                                Math.max(
                                    0,
                                    waveform.chirp_duration /
                                        waveform.chirp_period
                                )
                            )
                          : null;
                if (dutyCycle === null) {
                    return;
                }
                const txId = c.type === 'monostatic' ? c.txId : c.id;
                dutyCycleByTx.set(txId, dutyCycle);
            });
        });

        return dutyCycleByTx;
    }, [platforms, waveforms]);

    // Reset throttle and accumulated ranges when the scenario structure changes
    useEffect(() => {
        lastFetchTimeRef.current = 0;
        linkRangesRef.current.clear();
    }, [componentToPlatform]);

    // REPLACED: useFrame loop handles fetching instead of useEffect.
    // This prevents the cleanup-cancellation race condition where rapid
    // timeline updates would cancel the fetch request before it returned.
    useFrame(() => {
        // Access store state imperatively to avoid dependency staleness
        const state = useScenarioStore.getState();
        const simState = useSimulationProgressStore.getState();

        // 1. Concurrency & Sync Checks
        if (
            simState.isSimulating ||
            simState.isGeneratingKml ||
            state.isBackendSyncing ||
            isFetchingRef.current
        ) {
            return;
        }

        const now = Date.now();

        // 2. Throttle Check
        if (
            lastFetchTimeRef.current > 0 &&
            now - lastFetchTimeRef.current < FETCH_INTERVAL_MS
        ) {
            return;
        }

        // 3. Perform Fetch
        const fetchLinks = async () => {
            try {
                isFetchingRef.current = true;
                lastFetchTimeRef.current = Date.now();

                // Use the imperative time from the store for the most recent value
                const links = await invoke<RustVisualLink[]>(
                    'get_preview_links',
                    {
                        time: state.currentTime,
                    }
                );

                // Only update state if component is still mounted
                if (isMountedRef.current) {
                    const reconstructedLinks: LinkMetadata[] = links.map(
                        (l) => ({
                            link_type: TYPE_MAP[l.link_type],
                            quality: QUALITY_MAP[l.quality],
                            label: l.label,
                            source_id: l.source_id,
                            dest_id: l.dest_id,
                            origin_id: l.origin_id,
                            rcs: l.rcs,
                            actual_power_dbm: l.actual_power_dbm,
                            display_value: l.display_value,
                        })
                    );
                    setLinkMetadata(reconstructedLinks);

                    // Accumulate observed min/max per link (ranges only ever expand)
                    reconstructedLinks.forEach((m) => {
                        const key = `${m.link_type}_${m.source_id}_${m.dest_id}_${m.origin_id}`;
                        const existing = linkRangesRef.current.get(key);
                        if (existing) {
                            if (m.actual_power_dbm > -990) {
                                existing.powerMin = Math.min(
                                    existing.powerMin,
                                    m.actual_power_dbm
                                );
                                existing.powerMax = Math.max(
                                    existing.powerMax,
                                    m.actual_power_dbm
                                );
                            }
                            if (m.rcs >= 0) {
                                existing.rcsMin = Math.min(
                                    existing.rcsMin,
                                    m.rcs
                                );
                                existing.rcsMax = Math.max(
                                    existing.rcsMax,
                                    m.rcs
                                );
                            }
                        } else {
                            linkRangesRef.current.set(key, {
                                powerMin:
                                    m.actual_power_dbm > -990
                                        ? m.actual_power_dbm
                                        : Infinity,
                                powerMax:
                                    m.actual_power_dbm > -990
                                        ? m.actual_power_dbm
                                        : -Infinity,
                                rcsMin: m.rcs >= 0 ? m.rcs : Infinity,
                                rcsMax: m.rcs >= 0 ? m.rcs : -Infinity,
                            });
                        }
                    });
                }
            } catch (e) {
                console.error('Link preview error:', e);
            } finally {
                // Release lock
                if (isMountedRef.current) {
                    isFetchingRef.current = false;
                }
            }
        };

        void fetchLinks();
    });

    // Refs to each rendered LabelCluster's inner div for direct DOM transform updates
    const clusterDivRefs = useRef<(HTMLDivElement | null)[]>([]);
    // Latest clusters held in a ref so useFrame always reads the current value
    const clustersRef = useRef<
        Array<{ position: THREE.Vector3; links: RenderableLink[] }>
    >([]);
    // Scratch vector reused each frame to avoid per-frame allocation
    const ndcScratch = useMemo(() => new THREE.Vector3(), []);

    // Screen-space label placement: runs every frame, nudges overlapping cluster
    // divs apart via CSS transform so all labels remain visible.
    useFrame(() => {
        const clusterCount = clustersRef.current.length;
        if (clusterCount === 0) return;

        const w = gl.domElement.clientWidth;
        const h = gl.domElement.clientHeight;

        // Project all cluster anchor positions to screen space
        const sx = new Float32Array(clusterCount);
        const sy = new Float32Array(clusterCount);
        for (let i = 0; i < clusterCount; i++) {
            ndcScratch.copy(clustersRef.current[i].position).project(camera);
            sx[i] = (ndcScratch.x * 0.5 + 0.5) * w;
            sy[i] = (-ndcScratch.y * 0.5 + 0.5) * h;
        }

        // Iterative force-push: separate overlapping label centers
        const ox = new Float32Array(clusterCount); // x offsets (pixels)
        const oy = new Float32Array(clusterCount); // y offsets (pixels)
        for (let iter = 0; iter < 8; iter++) {
            for (let i = 0; i < clusterCount; i++) {
                for (let j = i + 1; j < clusterCount; j++) {
                    const dx = sx[j] + ox[j] - (sx[i] + ox[i]);
                    const dy = sy[j] + oy[j] - (sy[i] + oy[i]);
                    const distSq = dx * dx + dy * dy;
                    if (distSq >= LABEL_OVERLAP_PX * LABEL_OVERLAP_PX) continue;

                    const dist = Math.sqrt(distSq);
                    const push = (LABEL_OVERLAP_PX - dist) * 0.5;
                    // Push apart along the separation vector; use (1, 0) as fallback
                    const nx = dist > 0.5 ? dx / dist : 1;
                    const ny = dist > 0.5 ? dy / dist : 0;
                    ox[i] -= nx * push;
                    oy[i] -= ny * push;
                    ox[j] += nx * push;
                    oy[j] += ny * push;
                }
            }
        }

        // Apply computed offsets as CSS transforms on each cluster div
        for (let i = 0; i < clusterCount; i++) {
            const div = clusterDivRefs.current[i];
            if (!div) continue;
            div.style.transform = `translate(${ox[i].toFixed(1)}px, ${oy[i].toFixed(1)}px)`;
        }
    });

    // 2. High-Frequency Geometry Calculation (Runs every render/frame)
    const { clusters, flatLinks } = useMemo(() => {
        const calculatedLinks: RenderableLink[] = [];
        const clusterMap = new Map<
            string,
            { position: THREE.Vector3; links: RenderableLink[] }
        >();

        linkMetadata.forEach((meta) => {
            // Filter by type
            if (meta.link_type === 'monostatic' && !showLinkMonostatic) return;
            if (meta.link_type === 'illuminator' && !showLinkIlluminator)
                return;
            if (meta.link_type === 'scattered' && !showLinkScattered) return;
            if (meta.link_type === 'direct' && !showLinkDirect) return;

            const sourcePlat = componentToPlatform.get(meta.source_id);
            const destPlat = componentToPlatform.get(meta.dest_id);

            if (sourcePlat && destPlat) {
                const startPos = calculateInterpolatedPosition(
                    sourcePlat,
                    currentTime
                );
                const endPos = calculateInterpolatedPosition(
                    destPlat,
                    currentTime
                );
                const dist = startPos.distanceTo(endPos);

                const rangeKey = `${meta.link_type}_${meta.source_id}_${meta.dest_id}_${meta.origin_id}`;
                const range = linkRangesRef.current.get(rangeKey) ?? null;

                const renderLink: RenderableLink = {
                    ...meta,
                    start: startPos,
                    end: endPos,
                    distance: dist,
                    source_name:
                        componentToName.get(meta.source_id) ?? meta.source_id,
                    dest_name:
                        componentToName.get(meta.dest_id) ?? meta.dest_id,
                    origin_name:
                        componentToName.get(meta.origin_id) ?? meta.origin_id,
                    fmcwDutyCycle:
                        fmcwDutyCycleByTransmitter.get(meta.origin_id) ?? null,
                    range,
                };

                calculatedLinks.push(renderLink);

                // Clustering logic for labels - only if labels are enabled
                if (showLinkLabels) {
                    const mid = startPos.clone().lerp(endPos, 0.5);
                    // Round keys to cluster nearby labels
                    const key = `${mid.x.toFixed(1)}_${mid.y.toFixed(1)}_${mid.z.toFixed(1)}`;

                    if (!clusterMap.has(key)) {
                        clusterMap.set(key, { position: mid, links: [] });
                    }
                    clusterMap.get(key)!.links.push(renderLink);
                }
            }
        });

        const computedClusters = Array.from(clusterMap.values());
        clustersRef.current = computedClusters;
        return {
            clusters: computedClusters,
            flatLinks: calculatedLinks,
        };
    }, [
        linkMetadata,
        currentTime,
        componentToPlatform,
        componentToName,
        fmcwDutyCycleByTransmitter,
        showLinkLabels,
        showLinkMonostatic,
        showLinkIlluminator,
        showLinkScattered,
        showLinkDirect,
    ]);

    return (
        <group>
            {/* Render all lines individually for geometry */}
            {flatLinks.map((link, i) => (
                <LinkLine key={`line-${i}`} link={link} />
            ))}

            {/* Render aggregated label clusters */}
            {clusters.map((cluster, i) => (
                <LabelCluster
                    key={`cluster-${i}`}
                    position={cluster.position}
                    links={cluster.links}
                    divRef={(el) => {
                        clusterDivRefs.current[i] = el;
                    }}
                />
            ))}
        </group>
    );
}
