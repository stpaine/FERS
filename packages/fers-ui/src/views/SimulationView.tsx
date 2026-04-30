// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import KeyboardArrowDownIcon from '@mui/icons-material/KeyboardArrowDown';
import KeyboardArrowRightIcon from '@mui/icons-material/KeyboardArrowRight';
import MapIcon from '@mui/icons-material/Map';
import PlayCircleOutlineIcon from '@mui/icons-material/PlayCircleOutline';
import {
    Box,
    Button,
    Card,
    CardActions,
    CardContent,
    CircularProgress,
    Collapse,
    Grid,
    IconButton,
    LinearProgress,
    List,
    ListItem,
    ListItemText,
    Paper,
    Table,
    TableBody,
    TableCell,
    TableContainer,
    TableHead,
    TableRow,
    TextField,
    Typography,
} from '@mui/material';
import { invoke } from '@tauri-apps/api/core';
import { listen } from '@tauri-apps/api/event';
import { dirname, join } from '@tauri-apps/api/path';
import { open, save } from '@tauri-apps/plugin-dialog';
import React, { useEffect, useRef, useState } from 'react';
import { useScenarioStore } from '@/stores/scenarioStore';
import { getBlockingFmcwValidationMessage } from '@/stores/scenarioStore/fmcwValidation';
import {
    normalizeSimulationOutputMetadata,
    type RawSimulationOutputMetadata,
    type SimulationOutputFileMetadata,
    type SimulationProgressState,
    useSimulationProgressStore,
} from '@/stores/simulationProgressStore';
import {
    addSimulationProgressEvent,
    getSimulationProgressPercent,
    normalizeCompletedProgressSnapshot,
} from './simulationProgress';

export const SimulationView = React.memo(function SimulationView() {
    const [metadataExportPath, setMetadataExportPath] = useState<string | null>(
        null
    );
    const [expandedMetadataPaths, setExpandedMetadataPaths] = useState<
        Set<string>
    >(() => new Set());
    const isSimulating = useSimulationProgressStore(
        (state) => state.isSimulating
    );
    const isGeneratingKml = useSimulationProgressStore(
        (state) => state.isGeneratingKml
    );
    const setIsGeneratingKml = useSimulationProgressStore(
        (state) => state.setIsGeneratingKml
    );
    const simulationProgress = useSimulationProgressStore(
        (state) => state.simulationProgress
    );
    const simulationRunStatus = useSimulationProgressStore(
        (state) => state.simulationRunStatus
    );
    const simulationRunError = useSimulationProgressStore(
        (state) => state.simulationRunError
    );
    const simulationOutputMetadata = useSimulationProgressStore(
        (state) => state.simulationOutputMetadata
    );
    const startSimulationRun = useSimulationProgressStore(
        (state) => state.startSimulationRun
    );
    const setSimulationProgressSnapshot = useSimulationProgressStore(
        (state) => state.setSimulationProgressSnapshot
    );
    const setSimulationOutputMetadata = useSimulationProgressStore(
        (state) => state.setSimulationOutputMetadata
    );
    const completeSimulationRun = useSimulationProgressStore(
        (state) => state.completeSimulationRun
    );
    const failSimulationRun = useSimulationProgressStore(
        (state) => state.failSimulationRun
    );
    const showError = useScenarioStore((state) => state.showError);
    const showSuccess = useScenarioStore((state) => state.showSuccess);
    const scenarioFilePath = useScenarioStore(
        (state) => state.scenarioFilePath
    );
    const outputDirectory = useScenarioStore((state) => state.outputDirectory);
    const setOutputDirectory = useScenarioStore(
        (state) => state.setOutputDirectory
    );

    // Use a Ref to store incoming data to avoid triggering re-renders on every event
    const progressRef = useRef<Record<string, SimulationProgressState>>({});

    useEffect(() => {
        let animationFrameId: number | undefined;

        const flushProgress = () => {
            setSimulationProgressSnapshot({ ...progressRef.current });
        };

        // The update loop synchronizes the Ref data to the State at screen refresh rate
        const updateLoop = () => {
            if (useSimulationProgressStore.getState().isSimulating) {
                flushProgress();
                animationFrameId = requestAnimationFrame(updateLoop);
            }
        };

        const unlistenSimComplete = listen<void>('simulation-complete', () => {
            console.log('Simulation completed successfully.');
            progressRef.current = normalizeCompletedProgressSnapshot(
                progressRef.current
            );
            flushProgress();
            completeSimulationRun();
            if (animationFrameId !== undefined) {
                cancelAnimationFrame(animationFrameId);
            }
        });

        const unlistenSimError = listen<string>('simulation-error', (event) => {
            const errorMessage = `Simulation failed: ${event.payload}`;
            console.error(errorMessage);
            showError(errorMessage);
            flushProgress();
            failSimulationRun(errorMessage);
            if (animationFrameId !== undefined) {
                cancelAnimationFrame(animationFrameId);
            }
        });

        const unlistenSimProgress = listen<SimulationProgressState>(
            'simulation-progress',
            (event) => {
                progressRef.current = addSimulationProgressEvent(
                    progressRef.current,
                    event.payload
                );
            }
        );

        const unlistenOutputMetadata = listen<string>(
            'simulation-output-metadata',
            (event) => {
                try {
                    setSimulationOutputMetadata(
                        normalizeSimulationOutputMetadata(
                            JSON.parse(
                                event.payload
                            ) as RawSimulationOutputMetadata
                        )
                    );
                } catch (err) {
                    const errorMessage =
                        err instanceof Error ? err.message : String(err);
                    showError(
                        `Failed to decode simulation metadata: ${errorMessage}`
                    );
                }
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
            if (animationFrameId !== undefined) {
                cancelAnimationFrame(animationFrameId);
            }
            Promise.all([
                unlistenSimComplete,
                unlistenSimError,
                unlistenSimProgress,
                unlistenOutputMetadata,
                unlistenKmlComplete,
                unlistenKmlError,
            ]).then((unlisteners) => {
                unlisteners.forEach((unlisten) => unlisten());
            });
        };
    }, [
        isSimulating,
        setSimulationProgressSnapshot,
        setSimulationOutputMetadata,
        completeSimulationRun,
        failSimulationRun,
        setIsGeneratingKml,
        showError,
    ]);

    const getEffectiveOutputDir = async () => {
        if (outputDirectory) return outputDirectory;
        if (scenarioFilePath) {
            try {
                return await dirname(scenarioFilePath);
            } catch (e) {
                console.warn('Failed to get dirname of scenario file', e);
            }
        }
        return '.';
    };

    const handleSelectOutputDir = async () => {
        try {
            const selected = await open({
                directory: true,
                multiple: false,
                defaultPath: await getEffectiveOutputDir(),
            });
            if (typeof selected === 'string') {
                setOutputDirectory(selected);
            }
        } catch (err) {
            console.error('Failed to open directory dialog:', err);
        }
    };

    const handleRunSimulation = async () => {
        const scenarioState = useScenarioStore.getState();
        const validationMessage =
            getBlockingFmcwValidationMessage(scenarioState);
        if (validationMessage) {
            showError(`FMCW validation failed: ${validationMessage}`);
            return;
        }

        progressRef.current = {};
        setMetadataExportPath(null);
        startSimulationRun();
        try {
            // Ensure the C++ backend has the latest scenario from the UI
            const effectiveDir = await getEffectiveOutputDir();
            await invoke('set_output_directory', { dir: effectiveDir });

            await useScenarioStore.getState().syncBackend();
            await invoke('run_simulation');
        } catch (err) {
            const errorMessage =
                err instanceof Error ? err.message : String(err);
            console.error('Failed to invoke simulation:', errorMessage);
            showError(`Failed to start simulation: ${errorMessage}`);
            failSimulationRun(`Failed to start simulation: ${errorMessage}`);
        }
    };

    const handleGenerateKml = async () => {
        try {
            const scenarioState = useScenarioStore.getState();
            const validationMessage =
                getBlockingFmcwValidationMessage(scenarioState);
            if (validationMessage) {
                showError(`FMCW validation failed: ${validationMessage}`);
                return;
            }

            const effectiveDir = await getEffectiveOutputDir();

            // 1. Get the simulation name from the store
            const simName =
                scenarioState.globalParameters.simulation_name || 'scenario';

            // 2. Sanitize the name and append extension
            const suggestedFileName = `${simName.replace(/[^a-z0-9]/gi, '_')}.kml`;

            // 3. Join the directory and filename to create the pre-fill path
            const defaultPath = await join(effectiveDir, suggestedFileName);

            await invoke('set_output_directory', { dir: effectiveDir });

            const outputPath = await save({
                title: 'Save KML File',
                defaultPath: defaultPath,
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

    const hasProgress = Object.keys(simulationProgress).length > 0;
    const progressPanelVisible =
        isSimulating || hasProgress || simulationRunStatus === 'failed';
    const mainProgress = simulationProgress['main'];
    const progressHeading = mainProgress
        ? mainProgress.message
        : simulationRunStatus === 'completed'
          ? 'Simulation complete'
          : simulationRunStatus === 'failed'
            ? 'Simulation failed'
            : 'Preparing simulation...';
    const otherProgresses = Object.entries(simulationProgress)
        .filter(([key]) => key !== 'main')
        .sort((a, b) => a[0].localeCompare(b[0]));
    const mainProgressPercent = mainProgress
        ? getSimulationProgressPercent(mainProgress)
        : null;
    const renderProgressDetails = (
        details: SimulationProgressState['details']
    ) => {
        if (!details || details.length === 0) {
            return null;
        }

        return (
            <Box component="ul" sx={{ m: 0, mt: 1, pl: 2 }}>
                {details.map((detail) => {
                    const detailPercent = getSimulationProgressPercent(detail);
                    const detailSuffix =
                        detailPercent !== null
                            ? ` (${Math.round(detailPercent)}%)`
                            : detail.current > 0
                              ? ` (Chunk ${detail.current})`
                              : '';

                    return (
                        <Typography
                            component="li"
                            variant="caption"
                            color="text.secondary"
                            key={detail.id}
                        >
                            {detail.message}
                            {detailSuffix}
                        </Typography>
                    );
                })}
            </Box>
        );
    };
    const exportMetadataJson = async () => {
        try {
            const outputPath = await invoke<string>(
                'export_output_metadata_json'
            );
            setMetadataExportPath(outputPath);
            showSuccess(`Metadata JSON saved to ${outputPath}`);
        } catch (err) {
            const errorMessage =
                err instanceof Error ? err.message : String(err);
            showError(`Failed to export metadata JSON: ${errorMessage}`);
        }
    };
    const formatSampleRange = (start: number, end: number) =>
        `[${start}, ${end})`;
    const formatPulseLength = (
        minSamples: number,
        maxSamples: number,
        uniform: boolean
    ) => {
        if (minSamples === 0 && maxSamples === 0) {
            return '0';
        }
        return uniform ? String(minSamples) : `${minSamples} - ${maxSamples}`;
    };
    const formatMetric = (value: number) =>
        value.toLocaleString(undefined, { maximumSignificantDigits: 6 });
    const formatMetadataSamplingRates = () => {
        if (!simulationOutputMetadata) {
            return '';
        }
        if (typeof simulationOutputMetadata.sampling_rate === 'number') {
            return `${formatMetric(
                simulationOutputMetadata.sampling_rate
            )} samples/s`;
        }
        return `${simulationOutputMetadata.sampling_rates?.length ?? 0} sample rates`;
    };
    const formatFmcwMetadata = (
        fmcw: NonNullable<SimulationOutputFileMetadata['fmcw']>
    ) => {
        if (fmcw.waveform_shape === 'triangle') {
            return `triangle, B=${formatMetric(
                fmcw.chirp_bandwidth
            )}, T_c=${formatMetric(
                fmcw.chirp_duration
            )}, T_tri=${formatMetric(fmcw.triangle_period ?? 0)}`;
        }

        return `${fmcw.chirp_direction ?? 'up'}, B=${formatMetric(
            fmcw.chirp_bandwidth
        )}, T_c=${formatMetric(
            fmcw.chirp_duration
        )}, T_rep=${formatMetric(fmcw.chirp_period ?? 0)}`;
    };
    const formatPulseOrSegmentSummary = (
        file: SimulationOutputFileMetadata
    ) => {
        if (file.mode === 'pulsed') {
            return `${file.pulse_count} pulses, ${formatPulseLength(
                file.min_pulse_length_samples,
                file.max_pulse_length_samples,
                file.uniform_pulse_length
            )} samples`;
        }

        const segmentSummary = `${
            file.streaming_segments.length
        } streaming segments`;
        if (file.mode !== 'fmcw' || !file.fmcw) {
            if (file.mode === 'fmcw' && file.fmcw_sources.length > 0) {
                return `${segmentSummary}, ${file.fmcw_sources.length} FMCW sources`;
            }
            return segmentSummary;
        }

        return `${segmentSummary}, ${formatFmcwMetadata(file.fmcw)}`;
    };

    const toggleMetadataRow = (path: string) => {
        setExpandedMetadataPaths((current) => {
            const next = new Set(current);
            if (next.has(path)) {
                next.delete(path);
            } else {
                next.add(path);
            }
            return next;
        });
    };

    const renderFmcwDetails = (file: SimulationOutputFileMetadata) => {
        if (file.mode !== 'fmcw') {
            return null;
        }

        return (
            <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
                <Typography variant="body2">
                    File sample rate: {formatMetric(file.sampling_rate)}{' '}
                    samples/s
                </Typography>
                {file.fmcw_dechirp_mode &&
                    file.fmcw_dechirp_mode !== 'none' && (
                        <Typography variant="body2">
                            Dechirp: {file.fmcw_dechirp_mode},{' '}
                            {file.fmcw_dechirp_reference_source ?? 'none'}
                        </Typography>
                    )}
                {file.fmcw_sources.length === 0 ? (
                    <Typography variant="body2" color="text.secondary">
                        No FMCW source metadata was emitted for this receiver.
                    </Typography>
                ) : (
                    file.fmcw_sources.map((source) => (
                        <Box
                            key={`${source.transmitter_id}:${source.waveform_id}`}
                            sx={{ pl: 2 }}
                        >
                            <Typography variant="body2">
                                {source.transmitter_name} /{' '}
                                {source.waveform_name}:{' '}
                                {formatFmcwMetadata(source)}, f0=
                                {formatMetric(source.start_frequency_offset)}{' '}
                                Hz, rate={formatMetric(source.chirp_rate)}
                                Hz/s
                                {source.chirp_rate_signed !== undefined
                                    ? `, signed=${formatMetric(
                                          source.chirp_rate_signed
                                      )} Hz/s`
                                    : ''}
                                {source.chirp_count !== undefined
                                    ? `, chirps=${source.chirp_count}`
                                    : ''}
                                {source.triangle_count !== undefined
                                    ? `, triangles=${source.triangle_count}`
                                    : ''}
                            </Typography>
                            {source.segments.map((segment) => (
                                <Typography
                                    key={`${segment.start_time}:${segment.end_time}`}
                                    variant="caption"
                                    color="text.secondary"
                                    sx={{ display: 'block' }}
                                >
                                    [{formatMetric(segment.start_time)},{' '}
                                    {formatMetric(segment.end_time)}]:{' '}
                                    {segment.emitted_chirp_count !== undefined
                                        ? `${segment.emitted_chirp_count} chirps`
                                        : ''}
                                    {segment.emitted_triangle_count !==
                                    undefined
                                        ? `${segment.emitted_triangle_count} triangles`
                                        : ''}
                                </Typography>
                            ))}
                        </Box>
                    ))
                )}
            </Box>
        );
    };

    return (
        <Box
            sx={{ p: 4, height: '100%', overflowY: 'auto', contain: 'content' }}
        >
            <Typography variant="h4" gutterBottom>
                Simulation Runner
            </Typography>
            <Typography variant="body1" color="text.secondary" sx={{ mb: 4 }}>
                Execute the configured scenario or generate a geographical
                visualization. Ensure your scenario is fully configured before
                proceeding.
            </Typography>
            <Card elevation={0} sx={{ mb: 4 }}>
                <CardContent>
                    <Typography variant="h6" gutterBottom>
                        Output Settings
                    </Typography>
                    <Typography
                        variant="body2"
                        color="text.secondary"
                        sx={{ mb: 2 }}
                    >
                        Simulation results (.h5 files) and default KML exports
                        will be saved here.
                    </Typography>
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
                        <TextField
                            label="Output Directory"
                            variant="outlined"
                            size="small"
                            fullWidth
                            value={
                                outputDirectory ||
                                (scenarioFilePath
                                    ? 'Default (Scenario Directory)'
                                    : 'Default (Current Directory)')
                            }
                            slotProps={{
                                input: {
                                    readOnly: true,
                                },
                            }}
                        />
                        <Button
                            variant="outlined"
                            onClick={handleSelectOutputDir}
                            sx={{ whiteSpace: 'nowrap' }}
                        >
                            Browse...
                        </Button>
                        {outputDirectory && (
                            <Button
                                variant="text"
                                color="error"
                                onClick={() => setOutputDirectory(null)}
                            >
                                Reset
                            </Button>
                        )}
                    </Box>
                </CardContent>
            </Card>

            <Grid container spacing={4} sx={{ width: '100%' }}>
                {/* ... existing Grid items for Run Simulation and Generate KML ... */}
                <Grid size={{ xs: 12, md: 6 }}>
                    <Card elevation={0} sx={{ height: '100%' }}>
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
                    <Card elevation={0} sx={{ height: '100%' }}>
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

            {progressPanelVisible && (
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
                        {progressHeading}
                    </Typography>
                    {simulationRunStatus === 'failed' && simulationRunError && (
                        <Typography
                            variant="body2"
                            color="error"
                            sx={{ mb: 2, textAlign: 'center' }}
                        >
                            {simulationRunError}
                        </Typography>
                    )}
                    {mainProgressPercent !== null && (
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
                                    value={mainProgressPercent}
                                />
                            </Box>
                            <Box sx={{ minWidth: 40 }}>
                                <Typography
                                    variant="body2"
                                    color="text.secondary"
                                >{`${Math.round(mainProgressPercent)}%`}</Typography>
                            </Box>
                        </Box>
                    )}
                    {renderProgressDetails(mainProgress?.details)}

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
                                {otherProgresses.map(([key, prog]) => {
                                    const progressPercent =
                                        getSimulationProgressPercent(prog);
                                    const showChunkLabel =
                                        progressPercent === null &&
                                        prog.current > 0;

                                    return (
                                        <ListItem key={key}>
                                            <Box sx={{ width: '100%' }}>
                                                <Box
                                                    sx={{
                                                        display: 'flex',
                                                        alignItems: 'center',
                                                    }}
                                                >
                                                    <ListItemText
                                                        primary={prog.message}
                                                    />
                                                    {showChunkLabel && (
                                                        <Box
                                                            sx={{
                                                                width: '20%',
                                                                ml: 2,
                                                            }}
                                                        >
                                                            <Typography
                                                                variant="caption"
                                                                color="text.secondary"
                                                            >
                                                                Chunk{' '}
                                                                {prog.current}
                                                            </Typography>
                                                        </Box>
                                                    )}
                                                    {progressPercent !==
                                                        null && (
                                                        <Box
                                                            sx={{
                                                                width: '30%',
                                                                ml: 2,
                                                            }}
                                                        >
                                                            <LinearProgress
                                                                variant="determinate"
                                                                value={
                                                                    progressPercent
                                                                }
                                                            />
                                                        </Box>
                                                    )}
                                                </Box>
                                                {renderProgressDetails(
                                                    prog.details
                                                )}
                                            </Box>
                                        </ListItem>
                                    );
                                })}
                            </List>
                        </Box>
                    )}
                </Box>
            )}

            {simulationOutputMetadata && (
                <Card elevation={0} sx={{ mt: 4 }}>
                    <CardContent>
                        <Box
                            sx={{
                                display: 'flex',
                                justifyContent: 'space-between',
                                alignItems: 'center',
                                gap: 2,
                                mb: 2,
                            }}
                        >
                            <Box>
                                <Typography variant="h6">
                                    Output Data Metadata
                                </Typography>
                                <Typography
                                    variant="body2"
                                    color="text.secondary"
                                >
                                    {simulationOutputMetadata.files.length} HDF5
                                    output file
                                    {simulationOutputMetadata.files.length === 1
                                        ? ''
                                        : 's'}{' '}
                                    at {formatMetadataSamplingRates()}.
                                </Typography>
                            </Box>
                            <Button
                                variant="outlined"
                                onClick={exportMetadataJson}
                            >
                                Export JSON
                            </Button>
                        </Box>
                        {metadataExportPath && (
                            <Typography
                                variant="body2"
                                color="text.secondary"
                                sx={{ mb: 2, overflowWrap: 'anywhere' }}
                            >
                                Metadata JSON saved to {metadataExportPath}
                            </Typography>
                        )}

                        {simulationOutputMetadata.files.length === 0 ? (
                            <Typography color="text.secondary">
                                No HDF5 output files were generated for this
                                run.
                            </Typography>
                        ) : (
                            <TableContainer
                                component={Paper}
                                variant="outlined"
                            >
                                <Table size="small">
                                    <TableHead>
                                        <TableRow>
                                            <TableCell />
                                            <TableCell>Receiver</TableCell>
                                            <TableCell>Mode</TableCell>
                                            <TableCell align="right">
                                                Rate
                                            </TableCell>
                                            <TableCell align="right">
                                                Samples
                                            </TableCell>
                                            <TableCell>Sample Range</TableCell>
                                            <TableCell>Pulse/Segment</TableCell>
                                            <TableCell>File</TableCell>
                                        </TableRow>
                                    </TableHead>
                                    <TableBody>
                                        {simulationOutputMetadata.files.map(
                                            (file) => {
                                                const isExpanded =
                                                    expandedMetadataPaths.has(
                                                        file.path
                                                    );
                                                const canExpand =
                                                    file.mode === 'fmcw';
                                                return (
                                                    <React.Fragment
                                                        key={file.path}
                                                    >
                                                        <TableRow>
                                                            <TableCell>
                                                                {canExpand && (
                                                                    <IconButton
                                                                        size="small"
                                                                        onClick={() =>
                                                                            toggleMetadataRow(
                                                                                file.path
                                                                            )
                                                                        }
                                                                    >
                                                                        {isExpanded ? (
                                                                            <KeyboardArrowDownIcon fontSize="small" />
                                                                        ) : (
                                                                            <KeyboardArrowRightIcon fontSize="small" />
                                                                        )}
                                                                    </IconButton>
                                                                )}
                                                            </TableCell>
                                                            <TableCell>
                                                                {
                                                                    file.receiver_name
                                                                }
                                                            </TableCell>
                                                            <TableCell>
                                                                {file.mode}
                                                            </TableCell>
                                                            <TableCell align="right">
                                                                {formatMetric(
                                                                    file.sampling_rate
                                                                )}
                                                            </TableCell>
                                                            <TableCell align="right">
                                                                {
                                                                    file.total_samples
                                                                }
                                                            </TableCell>
                                                            <TableCell>
                                                                {formatSampleRange(
                                                                    file.sample_start,
                                                                    file.sample_end_exclusive
                                                                )}
                                                            </TableCell>
                                                            <TableCell>
                                                                {formatPulseOrSegmentSummary(
                                                                    file
                                                                )}
                                                            </TableCell>
                                                            <TableCell
                                                                sx={{
                                                                    maxWidth: 360,
                                                                    overflowWrap:
                                                                        'anywhere',
                                                                }}
                                                            >
                                                                {file.path}
                                                            </TableCell>
                                                        </TableRow>
                                                        {canExpand && (
                                                            <TableRow>
                                                                <TableCell
                                                                    colSpan={8}
                                                                    sx={{
                                                                        py: 0,
                                                                    }}
                                                                >
                                                                    <Collapse
                                                                        in={
                                                                            isExpanded
                                                                        }
                                                                        timeout="auto"
                                                                        unmountOnExit
                                                                    >
                                                                        <Box
                                                                            sx={{
                                                                                py: 2,
                                                                            }}
                                                                        >
                                                                            {renderFmcwDetails(
                                                                                file
                                                                            )}
                                                                        </Box>
                                                                    </Collapse>
                                                                </TableCell>
                                                            </TableRow>
                                                        )}
                                                    </React.Fragment>
                                                );
                                            }
                                        )}
                                    </TableBody>
                                </Table>
                            </TableContainer>
                        )}
                    </CardContent>
                </Card>
            )}
        </Box>
    );
});
