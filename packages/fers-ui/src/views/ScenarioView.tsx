// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import ChevronLeftIcon from '@mui/icons-material/ChevronLeft';
import {
    Alert,
    Box,
    Button,
    IconButton,
    Paper,
    Tooltip,
    Typography,
} from '@mui/material';
import { Canvas } from '@react-three/fiber';
import React, { useEffect, useRef, useState } from 'react';
import {
    Group,
    type GroupImperativeHandle,
    Panel,
    Separator,
} from 'react-resizable-panels';
import { type MapControls as MapControlsImpl } from 'three-stdlib';
import PropertyInspector from '@/components/PropertyInspector';
import { ScaleManager } from '@/components/ScaleManager';
import SceneTree from '@/components/SceneTree';
import Timeline from '@/components/Timeline';
import ViewControls from '@/components/ViewControls';
import { ViewportErrorBoundary } from '@/components/ViewportErrorBoundary';
import WorldView from '@/components/WorldView';
import { useScenarioStore } from '@/stores/scenarioStore';
import { fersColors } from '@/theme';
import {
    getWebGLSupportReport,
    resetWebGLSupportReportCache,
    type WebGLSupportReport,
} from '@/utils/webglSupport';

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

type ViewportState =
    | {
          phase: 'checking';
      }
    | {
          phase: 'ready';
          report: WebGLSupportReport;
      }
    | {
          phase: 'unsupported';
          report: WebGLSupportReport;
      }
    | {
          phase: 'renderer-error';
          report: WebGLSupportReport | null;
          errorMessage: string;
      };

function formatProbeSummary(report: WebGLSupportReport): string {
    return report.probes
        .map((probe) => `${probe.mode}: ${probe.reason}`)
        .join('\n');
}

function ScenarioViewportFallback({
    title,
    summary,
    diagnostics,
    onRetry,
    severity = 'warning',
}: {
    title: string;
    summary: string;
    diagnostics?: string;
    onRetry: () => void;
    severity?: 'warning' | 'error';
}) {
    return (
        <Box
            sx={{
                height: '100%',
                width: '100%',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                backgroundColor: fersColors.background.canvas,
                p: 3,
            }}
        >
            <Paper
                sx={{
                    width: 'min(720px, 100%)',
                    p: 3,
                    backgroundColor: fersColors.background.paper,
                }}
            >
                <Typography variant="h6" gutterBottom>
                    {title}
                </Typography>
                <Alert severity={severity} sx={{ mb: 2 }}>
                    {summary}
                </Alert>
                <Typography variant="body2" color="text.secondary">
                    The rest of the FERS UI remains available, but the 3D
                    Scenario view cannot be started on the current
                    browser/runtime/GPU configuration.
                </Typography>
                {diagnostics && (
                    <Box
                        component="pre"
                        sx={{
                            mt: 2,
                            mb: 0,
                            p: 2,
                            borderRadius: 1,
                            backgroundColor: 'rgba(2, 4, 8, 0.7)',
                            color: 'text.secondary',
                            fontFamily: 'monospace',
                            fontSize: '0.75rem',
                            overflowX: 'auto',
                            whiteSpace: 'pre-wrap',
                        }}
                    >
                        {diagnostics}
                    </Box>
                )}
                <Box
                    sx={{ mt: 2, display: 'flex', justifyContent: 'flex-end' }}
                >
                    <Button variant="outlined" onClick={onRetry}>
                        Retry 3D Startup
                    </Button>
                </Box>
            </Paper>
        </Box>
    );
}

/**
 * ScenarioView is the primary workbench for building and visualizing 3D scenes.
 */
export const ScenarioView = React.memo(function ScenarioView() {
    const isSimulating = useScenarioStore((state) => state.isSimulating);
    const [isInspectorCollapsed, setIsInspectorCollapsed] = useState(false);
    const panelGroupRef = useRef<GroupImperativeHandle>(null);
    const [viewportState, setViewportState] = useState<ViewportState>({
        phase: 'checking',
    });
    const [viewportResetKey, setViewportResetKey] = useState(0);

    // 1. Lift refs to this level to bridge Logic (Canvas) and UI (DOM)
    const controlsRef = useRef<MapControlsImpl>(null);
    const scaleLabelRef = useRef<HTMLDivElement>(null);
    const scaleBarRef = useRef<HTMLDivElement>(null);

    useEffect(() => {
        let cancelled = false;

        setViewportState({ phase: 'checking' });

        void getWebGLSupportReport()
            .then((report) => {
                if (cancelled) {
                    return;
                }

                if (!report.rendererSupported) {
                    console.warn(
                        'Scenario view disabled because WebGL2 startup checks failed.',
                        report
                    );
                    setViewportState({
                        phase: 'unsupported',
                        report,
                    });
                    return;
                }

                setViewportState({
                    phase: 'ready',
                    report,
                });
            })
            .catch((error) => {
                if (cancelled) {
                    return;
                }

                const errorMessage =
                    error instanceof Error
                        ? error.message
                        : 'Unknown WebGL startup error.';

                console.warn(
                    'Scenario view disabled because WebGL startup checks threw an unexpected error.',
                    error
                );
                setViewportState({
                    phase: 'renderer-error',
                    report: null,
                    errorMessage,
                });
            });

        return () => {
            cancelled = true;
        };
    }, [viewportResetKey]);

    const handleExpandInspector = () => {
        // Restore panels to their default sizes: [SceneTree, Main, Inspector]
        panelGroupRef.current?.setLayout({
            'scene-tree': 25,
            'main-content': 50,
            'property-inspector': 25,
        });
    };

    const handleViewportRetry = () => {
        resetWebGLSupportReportCache();
        setViewportResetKey((current) => current + 1);
    };

    const handleViewportError = (error: Error) => {
        console.error('Scenario view renderer startup failed.', error);
        setViewportState((current) => ({
            phase: 'renderer-error',
            report:
                current.phase === 'ready' || current.phase === 'unsupported'
                    ? current.report
                    : null,
            errorMessage: error.message,
        }));
    };

    const fallbackDiagnostics =
        viewportState.phase === 'unsupported'
            ? `${viewportState.report.summary}\n${formatProbeSummary(
                  viewportState.report
              )}`
            : viewportState.phase === 'renderer-error'
              ? [
                    viewportState.report?.summary,
                    viewportState.report
                        ? formatProbeSummary(viewportState.report)
                        : null,
                    `renderer: ${viewportState.errorMessage}`,
                ]
                    .filter(Boolean)
                    .join('\n')
              : undefined;

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
                            {viewportState.phase === 'checking' ? (
                                <ScenarioViewportFallback
                                    title="Checking 3D Renderer"
                                    summary="Running WebGL startup checks before loading the Scenario view."
                                    onRetry={handleViewportRetry}
                                />
                            ) : viewportState.phase === 'unsupported' ? (
                                <ScenarioViewportFallback
                                    title="3D Scenario View Unavailable"
                                    summary={viewportState.report.summary}
                                    diagnostics={fallbackDiagnostics}
                                    onRetry={handleViewportRetry}
                                />
                            ) : viewportState.phase === 'renderer-error' ? (
                                <ScenarioViewportFallback
                                    title="3D Scenario View Failed To Start"
                                    summary="Renderer startup failed after the initial WebGL probe. The viewport has been disabled to keep the rest of the app usable."
                                    diagnostics={fallbackDiagnostics}
                                    onRetry={handleViewportRetry}
                                    severity="error"
                                />
                            ) : (
                                <>
                                    <ViewportErrorBoundary
                                        resetKey={viewportResetKey}
                                        onError={handleViewportError}
                                    >
                                        <Canvas
                                            key={viewportResetKey}
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
                                                args={[
                                                    fersColors.background
                                                        .canvas,
                                                ]}
                                            />
                                            {/* Pass controlsRef down to WorldView */}
                                            <WorldView
                                                controlsRef={controlsRef}
                                            />

                                            {/* Logic-only component for Scale Bar calculations */}
                                            <ScaleManager
                                                controlsRef={controlsRef}
                                                labelRef={scaleLabelRef}
                                                barRef={scaleBarRef}
                                            />
                                        </Canvas>
                                    </ViewportErrorBoundary>

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
                                            textShadow:
                                                '0px 1px 2px rgba(0,0,0,0.8)',
                                        }}
                                    >
                                        <div
                                            ref={scaleLabelRef}
                                            style={{
                                                color: fersColors.text.primary,
                                                fontSize: '11px',
                                                fontWeight: 500,
                                                marginBottom: '4px',
                                                fontFamily:
                                                    'Roboto, sans-serif',
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
                                </>
                            )}
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
