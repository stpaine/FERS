// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import React, { useRef, useState } from 'react';
import { Box, IconButton, Tooltip } from '@mui/material';
import ChevronLeftIcon from '@mui/icons-material/ChevronLeft';
import { Canvas } from '@react-three/fiber';
import {
    Panel,
    Group,
    Separator,
    type GroupImperativeHandle,
} from 'react-resizable-panels';
import WorldView from '@/components/WorldView';
import SceneTree from '@/components/SceneTree';
import PropertyInspector from '@/components/PropertyInspector';
import Timeline from '@/components/Timeline';
import ViewControls from '@/components/ViewControls';
import { type MapControls as MapControlsImpl } from 'three-stdlib';
import { ScaleManager } from '@/components/ScaleManager';
import { useScenarioStore } from '@/stores/scenarioStore';
import { fersColors } from '@/theme';

/**
 * A styled resize handle for the resizable panels.
 */
function ResizeHandle() {
    return (
        <Separator>
            <Box
                sx={{
                    width: '2px',
                    height: '100%',
                    backgroundColor: 'divider',
                    transition: 'background-color 0.2s ease-in-out',
                    '&[data-resize-handle-state="drag"]': {
                        backgroundColor: 'primary.main',
                    },
                    '&:hover': {
                        backgroundColor: 'primary.light',
                    },
                }}
            />
        </Separator>
    );
}

/**
 * ScenarioView is the primary workbench for building and visualizing 3D scenes.
 */
export const ScenarioView = React.memo(function ScenarioView() {
    const isSimulating = useScenarioStore((state) => state.isSimulating);
    const [isInspectorCollapsed, setIsInspectorCollapsed] = useState(false);
    const panelGroupRef = useRef<GroupImperativeHandle>(null);

    // 1. Lift refs to this level to bridge Logic (Canvas) and UI (DOM)
    const controlsRef = useRef<MapControlsImpl>(null);
    const scaleLabelRef = useRef<HTMLDivElement>(null);
    const scaleBarRef = useRef<HTMLDivElement>(null);

    const handleExpandInspector = () => {
        // Restore panels to their default sizes: [SceneTree, Main, Inspector]
        panelGroupRef.current?.setLayout({
            'scene-tree': 25,
            'main-content': 50,
            'property-inspector': 25,
        });
    };

    return (
        <Box
            sx={{
                height: '100%',
                width: '100%',
                overflow: 'hidden',
                position: 'relative', // Establish positioning context
                pointerEvents: isSimulating ? 'none' : 'auto',
                opacity: isSimulating ? 0.5 : 1,
                transition: 'opacity 0.3s ease-in-out',
                userSelect: 'none',
                WebkitUserSelect: 'none',
            }}
        >
            <Group
                orientation="horizontal"
                groupRef={panelGroupRef}
                defaultLayout={{
                    'scene-tree': 25,
                    'main-content': 50,
                    'property-inspector': 25,
                }}
                style={{ height: '100%', width: '100%' }}
            >
                <Panel id="scene-tree" defaultSize={25} minSize={20}>
                    <SceneTree />
                </Panel>

                <ResizeHandle />

                <Panel id="main-content" minSize={30}>
                    <Box
                        sx={{
                            height: '100%',
                            display: 'flex',
                            flexDirection: 'column',
                            minWidth: 0, // Allow flex item to shrink below content size
                            overflow: 'hidden', // Prevent overflow
                        }}
                    >
                        <Box
                            sx={{
                                flex: 1,
                                position: 'relative',
                                minHeight: 0, // Allow flex item to shrink
                                overflow: 'hidden',
                                userSelect: 'none',
                            }}
                        >
                            <Canvas
                                shadows
                                camera={{
                                    position: [100, 100, 100],
                                    fov: 25,
                                    near: 0.1,
                                    far: 1e8,
                                }}
                                gl={{
                                    logarithmicDepthBuffer: true,
                                }}
                            >
                                <color
                                    attach="background"
                                    args={[fersColors.background.canvas]}
                                />
                                {/* Pass controlsRef down to WorldView */}
                                <WorldView controlsRef={controlsRef} />

                                {/* Logic-only component for Scale Bar calculations */}
                                <ScaleManager
                                    controlsRef={controlsRef}
                                    labelRef={scaleLabelRef}
                                    barRef={scaleBarRef}
                                />
                            </Canvas>

                            {/* View Controls (Top-Left) */}
                            <Box
                                sx={{
                                    position: 'absolute',
                                    top: 16,
                                    left: 16,
                                    zIndex: 1000,
                                }}
                            >
                                <ViewControls />
                            </Box>

                            {/* Scale Bar Overlay (Bottom-Left) */}
                            <Box
                                sx={{
                                    position: 'absolute',
                                    bottom: 16,
                                    left: 16,
                                    zIndex: 1000,
                                    pointerEvents: 'none',
                                    display: 'flex',
                                    flexDirection: 'column',
                                    alignItems: 'center',
                                    textShadow: '0px 1px 2px rgba(0,0,0,0.8)',
                                }}
                            >
                                <div
                                    ref={scaleLabelRef}
                                    style={{
                                        color: fersColors.text.primary,
                                        fontSize: '11px',
                                        fontWeight: 500,
                                        marginBottom: '4px',
                                        fontFamily: 'Roboto, sans-serif',
                                    }}
                                >
                                    -- m
                                </div>
                                <div
                                    ref={scaleBarRef}
                                    style={{
                                        height: '6px',
                                        border: `1px solid ${fersColors.text.primary}`,
                                        borderTop: 'none',
                                        width: '100px',
                                    }}
                                />
                            </Box>
                        </Box>
                        <Box
                            sx={{
                                height: 100,
                                flexShrink: 0,
                                borderTop: 1,
                                borderColor: 'divider',
                            }}
                        >
                            <Timeline />
                        </Box>
                    </Box>
                </Panel>

                <ResizeHandle />

                <Panel
                    id="property-inspector"
                    collapsible
                    collapsedSize={0}
                    onResize={(size) =>
                        setIsInspectorCollapsed(size.asPercentage === 0)
                    }
                    defaultSize={25}
                    minSize={5}
                >
                    <PropertyInspector />
                </Panel>
            </Group>

            {isInspectorCollapsed && (
                <Tooltip title="Show Properties">
                    <IconButton
                        onClick={handleExpandInspector}
                        sx={{
                            position: 'absolute',
                            right: 0,
                            top: '50%',
                            transform: 'translateY(-50%)',
                            zIndex: 1,
                            bgcolor: 'background.paper',
                            border: 1,
                            borderColor: 'divider',
                            borderRight: 'none',
                            borderTopRightRadius: 0,
                            borderBottomRightRadius: 0,
                            '&:hover': {
                                bgcolor: 'action.hover',
                            },
                        }}
                    >
                        <ChevronLeftIcon />
                    </IconButton>
                </Tooltip>
            )}
        </Box>
    );
});
