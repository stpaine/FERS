// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { useState, useEffect, useMemo, useRef } from 'react';
import { invoke } from '@tauri-apps/api/core';
import { Line, Html } from '@react-three/drei';
import { useFrame } from '@react-three/fiber';
import * as THREE from 'three';
import {
    useScenarioStore,
    Platform,
    calculateInterpolatedPosition,
} from '@/stores/scenarioStore';
import { Tooltip, Box, Typography } from '@mui/material';
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
}

// Derived object used for actual rendering
interface RenderableLink extends LinkMetadata {
    start: THREE.Vector3;
    end: THREE.Vector3;
    distance: number;
}

// Define the shape coming from Rust
interface RustVisualLink {
    link_type: number;
    quality: number;
    label: string;
    source_id: string;
    dest_id: string;
    origin_id: string;
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
                        <b>Path Segment:</b> {link.source_id} → {link.dest_id}
                    </Typography>
                    {link.link_type === 'scattered' && (
                        <Typography variant="caption" display="block">
                            <b>Illuminator:</b> {link.origin_id}
                        </Typography>
                    )}
                    <Typography variant="caption" display="block">
                        <b>Type:</b> {link.link_type}
                    </Typography>
                    <Typography variant="caption" display="block">
                        <b>Distance:</b> {(link.distance / 1000).toFixed(2)} km
                    </Typography>
                    <Typography variant="caption" display="block">
                        <b>Value:</b> {link.label}
                    </Typography>
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
}: {
    position: THREE.Vector3;
    links: RenderableLink[];
}) {
    return (
        <Html position={position} center zIndexRange={[100, 0]}>
            <div
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

/**
 * Renders radio frequency links between platforms.
 *
 * Performance Design:
 * 1. Metadata Fetching: Throttled to ~10 FPS. Fetches connectivity status and text labels.
 * 2. Geometry Rendering: Runs at 60 FPS (driven by store.currentTime).
 *
 * This ensures the lines "stick" to platforms smoothly while avoiding FFI/Text updates
 * overload on every frame.
 */
export default function LinkVisualizer() {
    const currentTime = useScenarioStore((state) => state.currentTime);
    const platforms = useScenarioStore((state) => state.platforms);
    const visibility = useScenarioStore((state) => state.visibility);
    const {
        showLinkLabels,
        showLinkMonostatic,
        showLinkIlluminator,
        showLinkScattered,
        showLinkDirect,
    } = visibility;

    // Store only the metadata (connectivity/labels), not positions
    const [linkMetadata, setLinkMetadata] = useState<LinkMetadata[]>([]);

    // Throttle control
    const lastFetchTimeRef = useRef<number>(0);
    const isFetchingRef = useRef<boolean>(false);
    const isMountedRef = useRef<boolean>(true);
    // Updated to 16ms to target 60FPS updates during playback
    const FETCH_INTERVAL_MS = 16;

    // Track component mount status to prevent setting state on unmounted component
    useEffect(() => {
        isMountedRef.current = true;
        return () => {
            isMountedRef.current = false;
        };
    }, []);

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

    // Reset throttle when the scenario structure changes
    useEffect(() => {
        lastFetchTimeRef.current = 0;
    }, [componentToPlatform]);

    // REPLACED: useFrame loop handles fetching instead of useEffect.
    // This prevents the cleanup-cancellation race condition where rapid
    // timeline updates would cancel the fetch request before it returned.
    useFrame(() => {
        // Access store state imperatively to avoid dependency staleness
        const state = useScenarioStore.getState();

        // 1. Concurrency & Sync Checks
        if (state.isBackendSyncing || isFetchingRef.current) {
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
                        })
                    );
                    setLinkMetadata(reconstructedLinks);
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

                const renderLink: RenderableLink = {
                    ...meta,
                    start: startPos,
                    end: endPos,
                    distance: dist,
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

        return {
            clusters: Array.from(clusterMap.values()),
            flatLinks: calculatedLinks,
        };
    }, [
        linkMetadata,
        currentTime,
        componentToPlatform,
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
                />
            ))}
        </group>
    );
}
