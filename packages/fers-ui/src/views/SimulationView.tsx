// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import React, { useState, useEffect, useRef } from 'react';
import {
    Box,
    Typography,
    Grid,
    Card,
    CardContent,
    CardActions,
    Button,
    CircularProgress,
    Fade,
    LinearProgress,
    List,
    ListItem,
    ListItemText,
} from '@mui/material';
import PlayCircleOutlineIcon from '@mui/icons-material/PlayCircleOutline';
import MapIcon from '@mui/icons-material/Map';
import { useScenarioStore } from '@/stores/scenarioStore';
import { invoke } from '@tauri-apps/api/core';
import { save } from '@tauri-apps/plugin-dialog';
import { listen } from '@tauri-apps/api/event';

interface ProgressState {
    message: string;
    current: number;
    total: number;
}

export const SimulationView = React.memo(function SimulationView() {
    const isSimulating = useScenarioStore((state) => state.isSimulating);
    const setIsSimulating = useScenarioStore((state) => state.setIsSimulating);
    const showError = useScenarioStore((state) => state.showError);
    const [isGeneratingKml, setIsGeneratingKml] = useState(false);

    // Use a Ref to store incoming data to avoid triggering re-renders on every event
    const progressRef = useRef<Record<string, ProgressState>>({});
    // Local state to trigger actual re-renders throttled by rAF
    const [displayProgress, setDisplayProgress] = useState<
        Record<string, ProgressState>
    >({});

    useEffect(() => {
        let animationFrameId: number;

        // The update loop synchronizes the Ref data to the State at screen refresh rate
        const updateLoop = () => {
            if (useScenarioStore.getState().isSimulating) {
                setDisplayProgress({ ...progressRef.current });
                animationFrameId = requestAnimationFrame(updateLoop);
            }
        };

        const unlistenSimComplete = listen<void>('simulation-complete', () => {
            console.log('Simulation completed successfully.');
            setIsSimulating(false);
            progressRef.current = {};
            setDisplayProgress({});
            cancelAnimationFrame(animationFrameId);
        });

        const unlistenSimError = listen<string>('simulation-error', (event) => {
            const errorMessage = `Simulation failed: ${event.payload}`;
            console.error(errorMessage);
            showError(errorMessage);
            setIsSimulating(false);
            progressRef.current = {};
            setDisplayProgress({});
            cancelAnimationFrame(animationFrameId);
        });

        const unlistenSimProgress = listen<ProgressState>(
            'simulation-progress',
            (event) => {
                const { message } = event.payload;
                let key = message;

                // Grouping logic to keep the list clean
                if (
                    message.startsWith('Simulating') ||
                    message.startsWith('Initializing')
                ) {
                    key = 'main';
                } else if (
                    message.startsWith('Finalizing') ||
                    message.startsWith('Exporting') ||
                    message.startsWith('Rendering') ||
                    message.startsWith('Applying') ||
                    message.startsWith('Writing') ||
                    message.startsWith('Finished')
                ) {
                    // Extract component name from format "Action {Name}: ..." or "Action {Name}"
                    // e.g. "Exporting Receiver1: Chunk 50" -> key "Receiver1"
                    const parts = message.split(/[:\s]+/);
                    // Simple heuristic: 2nd word is often the name in our C++ format
                    if (parts.length >= 2) {
                        key = parts[1];
                    }
                }

                progressRef.current[key] = event.payload;
            }
        );

        const unlistenKmlComplete = listen<string>(
            'kml-generation-complete',
            (event) => {
                console.log('KML generated successfully at:', event.payload);
                setIsGeneratingKml(false);
            }
        );

        const unlistenKmlError = listen<string>(
            'kml-generation-error',
            (event) => {
                const errorMessage = `KML generation failed: ${event.payload}`;
                console.error(errorMessage);
                showError(errorMessage);
                setIsGeneratingKml(false);
            }
        );

        // Start the UI update loop if we are simulating
        if (isSimulating) {
            updateLoop();
        }

        return () => {
            cancelAnimationFrame(animationFrameId);
            Promise.all([
                unlistenSimComplete,
                unlistenSimError,
                unlistenSimProgress,
                unlistenKmlComplete,
                unlistenKmlError,
            ]).then((unlisteners) => {
                unlisteners.forEach((unlisten) => unlisten());
            });
        };
    }, [isSimulating, setIsSimulating, showError]);

    const handleRunSimulation = async () => {
        progressRef.current = {};
        setDisplayProgress({});
        setIsSimulating(true);
        try {
            // Ensure the C++ backend has the latest scenario from the UI
            await useScenarioStore.getState().syncBackend();
            await invoke('run_simulation');
        } catch (err) {
            const errorMessage =
                err instanceof Error ? err.message : String(err);
            console.error('Failed to invoke simulation:', errorMessage);
            showError(`Failed to start simulation: ${errorMessage}`);
            setIsSimulating(false); // Stop on invocation failure
        }
    };

    const handleGenerateKml = async () => {
        try {
            const outputPath = await save({
                title: 'Save KML File',
                filters: [{ name: 'KML File', extensions: ['kml'] }],
            });

            if (outputPath) {
                setIsGeneratingKml(true);
                // Ensure the C++ backend has the latest scenario from the UI
                await useScenarioStore.getState().syncBackend();
                await invoke('generate_kml', { outputPath });
            }
        } catch (err) {
            const errorMessage =
                err instanceof Error ? err.message : String(err);
            console.error('Failed to invoke KML generation:', errorMessage);
            showError(`Failed to start KML generation: ${errorMessage}`);
            setIsGeneratingKml(false); // Stop on invocation failure
        }
    };

    const mainProgress = displayProgress['main'];
    const otherProgresses = Object.entries(displayProgress)
        .filter(([key]) => key !== 'main')
        .sort((a, b) => a[0].localeCompare(b[0]));

    return (
        <Box sx={{ p: 4, height: '100%', overflowY: 'auto' }}>
            <Typography variant="h4" gutterBottom>
                Simulation Runner
            </Typography>
            <Typography variant="body1" color="text.secondary" sx={{ mb: 4 }}>
                Execute the configured scenario or generate a geographical
                visualization. Ensure your scenario is fully configured before
                proceeding.
            </Typography>

            <Grid container spacing={4} sx={{ width: '100%' }}>
                <Grid size={{ xs: 12, md: 6 }}>
                    <Card sx={{ height: '100%' }}>
                        <CardContent>
                            <Typography variant="h5" component="div">
                                Run Full Simulation
                            </Typography>
                            <Typography sx={{ mt: 1.5 }} color="text.secondary">
                                Executes the entire simulation based on the
                                current scenario settings. This is a
                                computationally intensive process that will
                                generate output files.
                            </Typography>
                        </CardContent>
                        <CardActions sx={{ p: 2 }}>
                            <Button
                                variant="contained"
                                size="large"
                                startIcon={
                                    isSimulating ? (
                                        <CircularProgress
                                            size={24}
                                            color="inherit"
                                        />
                                    ) : (
                                        <PlayCircleOutlineIcon />
                                    )
                                }
                                disabled={isSimulating || isGeneratingKml}
                                onClick={handleRunSimulation}
                            >
                                {isSimulating ? 'Running...' : 'Run Simulation'}
                            </Button>
                        </CardActions>
                    </Card>
                </Grid>
                <Grid size={{ xs: 12, md: 6 }}>
                    <Card sx={{ height: '100%' }}>
                        <CardContent>
                            <Typography variant="h5" component="div">
                                Generate KML
                            </Typography>
                            <Typography sx={{ mt: 1.5 }} color="text.secondary">
                                Creates a KML file from the scenario&apos;s
                                platform motion paths and antenna pointings.
                                This allows for quick visualization in
                                applications like Google Earth without running
                                the full signal-level simulation.
                            </Typography>
                        </CardContent>
                        <CardActions sx={{ p: 2 }}>
                            <Button
                                variant="outlined"
                                size="large"
                                startIcon={
                                    isGeneratingKml ? (
                                        <CircularProgress
                                            size={24}
                                            color="inherit"
                                        />
                                    ) : (
                                        <MapIcon />
                                    )
                                }
                                disabled={isSimulating || isGeneratingKml}
                                onClick={handleGenerateKml}
                            >
                                {isGeneratingKml
                                    ? 'Generating...'
                                    : 'Generate KML'}
                            </Button>
                        </CardActions>
                    </Card>
                </Grid>
            </Grid>

            <Fade in={isSimulating}>
                <Box
                    sx={{
                        mt: 4,
                        p: 2,
                        backgroundColor: 'action.hover',
                        borderRadius: 1,
                    }}
                >
                    {/* Main Simulation Progress */}
                    <Typography
                        variant="h6"
                        sx={{ mb: 1, textAlign: 'center' }}
                    >
                        {mainProgress
                            ? mainProgress.message
                            : 'Preparing simulation...'}
                    </Typography>
                    {mainProgress && mainProgress.total > 0 && (
                        <Box
                            sx={{
                                display: 'flex',
                                alignItems: 'center',
                                mt: 2,
                                mb: 2,
                            }}
                        >
                            <Box sx={{ width: '100%', mr: 1 }}>
                                <LinearProgress
                                    variant="determinate"
                                    value={
                                        (mainProgress.current /
                                            mainProgress.total) *
                                        100
                                    }
                                />
                            </Box>
                            <Box sx={{ minWidth: 40 }}>
                                <Typography
                                    variant="body2"
                                    color="text.secondary"
                                >{`${Math.round(
                                    (mainProgress.current /
                                        mainProgress.total) *
                                        100
                                )}%`}</Typography>
                            </Box>
                        </Box>
                    )}

                    {/* Finalizer Threads List */}
                    {otherProgresses.length > 0 && (
                        <Box
                            sx={{
                                mt: 2,
                                borderTop: 1,
                                borderColor: 'divider',
                                pt: 2,
                            }}
                        >
                            <Typography
                                variant="subtitle2"
                                color="text.secondary"
                            >
                                Exporting Data:
                            </Typography>
                            <List dense>
                                {otherProgresses.map(([key, prog]) => (
                                    <ListItem key={key}>
                                        <ListItemText
                                            primary={prog.message}
                                            secondary={
                                                prog.total > 0 &&
                                                prog.total !== 100
                                                    ? 'Processing...'
                                                    : ''
                                            }
                                        />
                                        {prog.total > 0 &&
                                            prog.total !== 100 && (
                                                /* For chunks, we often don't know the exact total ahead of time,
                                               so a determinate bar might jump, but passing a large number
                                               or indeterminate looks better than flickering 100% */
                                                <Box
                                                    sx={{ width: '20%', ml: 2 }}
                                                >
                                                    <Typography
                                                        variant="caption"
                                                        color="text.secondary"
                                                    >
                                                        Chunk {prog.current}
                                                    </Typography>
                                                </Box>
                                            )}
                                        {prog.total === 100 && (
                                            <Box sx={{ width: '30%', ml: 2 }}>
                                                <LinearProgress
                                                    variant="determinate"
                                                    value={prog.current}
                                                />
                                            </Box>
                                        )}
                                    </ListItem>
                                ))}
                            </List>
                        </Box>
                    )}
                </Box>
            </Fade>
        </Box>
    );
});
