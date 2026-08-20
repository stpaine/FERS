// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import FilterListIcon from '@mui/icons-material/FilterList';
import PlayCircleOutlineIcon from '@mui/icons-material/PlayCircleOutline';
import SaveAltIcon from '@mui/icons-material/SaveAlt';
import StopCircleIcon from '@mui/icons-material/StopCircle';
import {
    Alert,
    Box,
    Button,
    Chip,
    FormControl,
    Grid,
    InputLabel,
    MenuItem,
    Paper,
    Select,
    Stack,
    Switch,
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
import { dirname } from '@tauri-apps/api/path';
import React, {
    useCallback,
    useEffect,
    useMemo,
    useRef,
    useState,
} from 'react';
import {
    totalLatePacketCount,
    Vita49LatePacketInfo,
} from '@/components/Vita49LatePacketInfo';
import { useScenarioStore } from '@/stores/scenarioStore';
import { getBlockingFmcwValidationMessage } from '@/stores/scenarioStore/fmcwValidation';
import {
    normalizeSimulationOutputMetadata,
    type RawSimulationOutputMetadata,
} from '@/stores/simulationProgressStore';
import {
    deriveExpectedVita49Streams,
    mergeVita49StreamRows,
    toVita49BackendConfig,
    useVita49StreamingStore,
    type Vita49PacketTraceEvent,
    type Vita49StreamStatsEvent,
    type Vita49TelemetryPoll,
    type Vita49Timestamp,
    validateVita49Config,
} from '@/stores/vita49StreamingStore';

const TELEMETRY_POLL_INTERVAL_MS = 100;
const TELEMETRY_POLL_ERROR_INTERVAL_MS = 250;
const TELEMETRY_POLL_PACKET_LIMIT = 1000;
const PACKET_TRACE_TABLE_HEIGHT = 420;
const PACKET_TRACE_ROW_HEIGHT = 36;
const PACKET_TRACE_OVERSCAN = 8;

const formatMetric = (value: number | null | undefined) =>
    value === null || value === undefined
        ? '-'
        : value.toLocaleString(undefined, { maximumSignificantDigits: 6 });

const formatExact = (value: string | number | null | undefined) =>
    value === null || value === undefined ? '-' : String(value);

const formatSeconds = (value: number | null | undefined) =>
    value === null || value === undefined
        ? '-'
        : `${value.toLocaleString(undefined, { maximumSignificantDigits: 9 })} s`;

const formatVita49Timestamp = (
    timestamp: Vita49Timestamp | null | undefined
) => {
    if (!timestamp) {
        return '-';
    }

    const second = Number(timestamp.integer_seconds);
    if (!Number.isFinite(second)) {
        return '-';
    }

    const base = new Date(second * 1000).toISOString().replace('.000Z', '');
    const fractionalPicoseconds = Math.trunc(timestamp.fractional_picoseconds)
        .toString()
        .padStart(12, '0')
        .replace(/0+$/, '');
    return fractionalPicoseconds
        ? `${base}.${fractionalPicoseconds}Z`
        : `${base}Z`;
};

const formatSimulationSpan = (
    start: number | null | undefined,
    end: number | null | undefined
) => {
    const first = formatSeconds(start);
    const last = formatSeconds(end);
    return first === '-' && last === '-' ? '-' : `${first} - ${last}`;
};

const formatTimestampSpan = (
    start: Vita49Timestamp | null | undefined,
    end: Vita49Timestamp | null | undefined
) => {
    const first = formatVita49Timestamp(start);
    const last = formatVita49Timestamp(end);
    return first === '-' && last === '-' ? '-' : `${first} - ${last}`;
};

const formatStreamId = (streamId: number | null | undefined) =>
    streamId === null || streamId === undefined
        ? '-'
        : `0x${streamId.toString(16).toUpperCase().padStart(8, '0')}`;

export const Vita49StreamingView = React.memo(function Vita49StreamingView() {
    const config = useVita49StreamingStore((state) => state.config);
    const runState = useVita49StreamingStore((state) => state.runState);
    const expectedStreams = useVita49StreamingStore(
        (state) => state.expectedStreams
    );
    const streamStats = useVita49StreamingStore((state) => state.streamStats);
    const packetTrace = useVita49StreamingStore((state) => state.packetTrace);
    const omittedPacketTraceEvents = useVita49StreamingStore(
        (state) => state.omittedPacketTraceEvents
    );
    const finalVita49Metadata = useVita49StreamingStore(
        (state) => state.finalVita49Metadata
    );
    const error = useVita49StreamingStore((state) => state.error);
    const setConfig = useVita49StreamingStore((state) => state.setConfig);
    const startRun = useVita49StreamingStore((state) => state.startRun);
    const markStopping = useVita49StreamingStore((state) => state.markStopping);
    const markDraining = useVita49StreamingStore((state) => state.markDraining);
    const setStreamStats = useVita49StreamingStore(
        (state) => state.setStreamStats
    );
    const appendPacketBatch = useVita49StreamingStore(
        (state) => state.appendPacketBatch
    );
    const completeRun = useVita49StreamingStore((state) => state.completeRun);
    const cancelRun = useVita49StreamingStore((state) => state.cancelRun);
    const failRun = useVita49StreamingStore((state) => state.failRun);
    const showError = useScenarioStore((state) => state.showError);
    const showWarning = useScenarioStore((state) => state.showWarning);
    const showSuccess = useScenarioStore((state) => state.showSuccess);
    const scenarioFilePath = useScenarioStore(
        (state) => state.scenarioFilePath
    );
    const outputDirectory = useScenarioStore((state) => state.outputDirectory);

    const [metadataExportPath, setMetadataExportPath] = useState<string | null>(
        null
    );
    const [streamFilter, setStreamFilter] = useState('all');
    const [packetKindFilter, setPacketKindFilter] = useState('all');
    const [droppedOnly, setDroppedOnly] = useState(false);
    const [overRangeOnly, setOverRangeOnly] = useState(false);
    const [sampleLossOnly, setSampleLossOnly] = useState(false);
    const [packetTraceScrollTop, setPacketTraceScrollTop] = useState(0);
    const packetTraceContainerRef = useRef<HTMLDivElement | null>(null);
    const pendingTelemetryRef = useRef<{
        stats: Vita49StreamStatsEvent | null;
        packets: Vita49PacketTraceEvent[];
        omittedPacketTraceEvents: number;
    }>({
        stats: null,
        packets: [],
        omittedPacketTraceEvents: 0,
    });
    const telemetryFlushFrameRef = useRef<number | null>(null);
    const isRunning =
        runState === 'running' ||
        runState === 'stopping' ||
        runState === 'draining';
    const configErrors = validateVita49Config(config);

    const appendTelemetry = useCallback(
        (telemetry: Vita49TelemetryPoll) => {
            if (telemetry.stats) {
                setStreamStats(telemetry.stats);
            }
            if (
                telemetry.packets.length > 0 ||
                telemetry.omitted_packet_trace_events > 0
            ) {
                appendPacketBatch(
                    telemetry.packets,
                    telemetry.omitted_packet_trace_events
                );
            }
        },
        [appendPacketBatch, setStreamStats]
    );

    const flushPendingTelemetry = useCallback(() => {
        if (telemetryFlushFrameRef.current !== null) {
            window.cancelAnimationFrame(telemetryFlushFrameRef.current);
            telemetryFlushFrameRef.current = null;
        }

        const pending = pendingTelemetryRef.current;
        pendingTelemetryRef.current = {
            stats: null,
            packets: [],
            omittedPacketTraceEvents: 0,
        };

        if (pending.stats) {
            setStreamStats(pending.stats);
        }
        if (
            pending.packets.length > 0 ||
            pending.omittedPacketTraceEvents > 0
        ) {
            appendPacketBatch(
                pending.packets,
                pending.omittedPacketTraceEvents
            );
        }
    }, [appendPacketBatch, setStreamStats]);

    const drainAvailableTelemetry = useCallback(async () => {
        flushPendingTelemetry();
        let hasMore = false;
        do {
            const telemetry = await invoke<Vita49TelemetryPoll>(
                'poll_vita49_telemetry',
                { maxPackets: TELEMETRY_POLL_PACKET_LIMIT }
            );
            appendTelemetry(telemetry);
            hasMore = telemetry.has_more;
        } while (hasMore);
        flushPendingTelemetry();
    }, [appendTelemetry, flushPendingTelemetry]);

    useEffect(() => {
        let active = true;
        const unlisteners = Promise.all([
            listen<string>('vita49-output-metadata', (event) => {
                if (!active) return;
                const metadata = normalizeSimulationOutputMetadata(
                    JSON.parse(event.payload) as RawSimulationOutputMetadata
                );
                void drainAvailableTelemetry().finally(() => {
                    if (active) completeRun(metadata);
                });
            }),
            listen<string>('vita49-stream-complete', (event) => {
                if (!active) return;
                const metadata = normalizeSimulationOutputMetadata(
                    JSON.parse(event.payload) as RawSimulationOutputMetadata
                );
                void drainAvailableTelemetry().finally(() => {
                    if (active) completeRun(metadata);
                });
            }),
            listen<string>('vita49-stream-draining', () => {
                if (!active) return;
                markDraining();
            }),
            listen<string>('vita49-stream-cancelled', (event) => {
                if (!active) return;
                const metadata = normalizeSimulationOutputMetadata(
                    JSON.parse(event.payload) as RawSimulationOutputMetadata
                );
                void drainAvailableTelemetry().finally(() => {
                    if (active) cancelRun(metadata);
                });
            }),
            listen<string>('vita49-stream-error', (event) => {
                if (!active) return;
                failRun(event.payload);
                showError(`VITA49 streaming failed: ${event.payload}`);
            }),
        ]);

        return () => {
            active = false;
            unlisteners.then((listeners) =>
                listeners.forEach((unlisten) => unlisten())
            );
        };
    }, [
        cancelRun,
        completeRun,
        drainAvailableTelemetry,
        failRun,
        markDraining,
        showError,
    ]);

    useEffect(() => {
        if (!isRunning) {
            return;
        }

        let cancelled = false;
        let pollTimer: ReturnType<typeof setTimeout> | null = null;

        const scheduleFlush = () => {
            if (telemetryFlushFrameRef.current !== null) {
                return;
            }
            telemetryFlushFrameRef.current = window.requestAnimationFrame(
                flushPendingTelemetry
            );
        };

        const queueTelemetry = (telemetry: Vita49TelemetryPoll) => {
            const pending = pendingTelemetryRef.current;
            pending.stats = telemetry.stats ?? pending.stats;
            pending.packets.push(...telemetry.packets);
            pending.omittedPacketTraceEvents +=
                telemetry.omitted_packet_trace_events;
            scheduleFlush();
        };

        const pollTelemetry = async () => {
            if (cancelled) {
                return;
            }

            try {
                let hasMore = false;
                do {
                    const telemetry = await invoke<Vita49TelemetryPoll>(
                        'poll_vita49_telemetry',
                        { maxPackets: TELEMETRY_POLL_PACKET_LIMIT }
                    );
                    queueTelemetry(telemetry);
                    hasMore = telemetry.has_more;
                } while (!cancelled && hasMore);

                if (!cancelled) {
                    pollTimer = setTimeout(
                        pollTelemetry,
                        TELEMETRY_POLL_INTERVAL_MS
                    );
                }
            } catch (err) {
                console.error('Failed to poll VITA49 telemetry:', err);
                if (!cancelled) {
                    pollTimer = setTimeout(
                        pollTelemetry,
                        TELEMETRY_POLL_ERROR_INTERVAL_MS
                    );
                }
            }
        };

        void pollTelemetry();

        return () => {
            cancelled = true;
            if (pollTimer !== null) {
                clearTimeout(pollTimer);
            }
            flushPendingTelemetry();
        };
    }, [flushPendingTelemetry, isRunning]);

    const streamRows = useMemo(
        () => mergeVita49StreamRows(expectedStreams, streamStats),
        [expectedStreams, streamStats]
    );

    const aggregate = useMemo(
        () =>
            streamRows.reduce(
                (acc, row) => ({
                    packets: acc.packets + row.packetsEmitted,
                    samples: acc.samples + row.samplesEmitted,
                    drops: acc.drops + row.packetsDropped,
                    lateData: acc.lateData + row.lateDataPacketCount,
                    lateContext: acc.lateContext + row.lateContextPacketCount,
                    overRange: acc.overRange + row.overRangeCount,
                    context: acc.context + row.contextPackets,
                }),
                {
                    packets: 0,
                    samples: 0,
                    drops: 0,
                    lateData: 0,
                    lateContext: 0,
                    overRange: 0,
                    context: 0,
                }
            ),
        [streamRows]
    );

    const streamIdOptions = useMemo(() => {
        const ids = new Set<number>();
        for (const row of streamRows) {
            if (row.streamId !== null) ids.add(row.streamId);
        }
        for (const packet of packetTrace) {
            if (packet.stream_id) ids.add(packet.stream_id);
        }
        return Array.from(ids).sort((a, b) => a - b);
    }, [packetTrace, streamRows]);

    const filteredPackets = useMemo(
        () =>
            packetTrace.filter((packet) => {
                if (
                    streamFilter !== 'all' &&
                    packet.stream_id !== Number(streamFilter)
                ) {
                    return false;
                }
                if (packetKindFilter === 'data' && !packet.data_packet) {
                    return false;
                }
                if (packetKindFilter === 'context' && !packet.context_packet) {
                    return false;
                }
                if (droppedOnly && !packet.dropped) return false;
                if (overRangeOnly && !packet.over_range) return false;
                if (sampleLossOnly && !packet.sample_loss) return false;
                return true;
            }),
        [
            droppedOnly,
            overRangeOnly,
            packetKindFilter,
            packetTrace,
            sampleLossOnly,
            streamFilter,
        ]
    );

    useEffect(() => {
        setPacketTraceScrollTop(0);
        if (packetTraceContainerRef.current) {
            packetTraceContainerRef.current.scrollTop = 0;
        }
    }, [
        droppedOnly,
        overRangeOnly,
        packetKindFilter,
        sampleLossOnly,
        streamFilter,
    ]);

    const packetTraceWindow = useMemo(() => {
        const rawStart = Math.max(
            0,
            Math.floor(packetTraceScrollTop / PACKET_TRACE_ROW_HEIGHT) -
                PACKET_TRACE_OVERSCAN
        );
        const start = Math.min(filteredPackets.length, rawStart);
        const end = Math.min(
            filteredPackets.length,
            Math.ceil(
                (packetTraceScrollTop + PACKET_TRACE_TABLE_HEIGHT) /
                    PACKET_TRACE_ROW_HEIGHT
            ) + PACKET_TRACE_OVERSCAN
        );
        return {
            start,
            end,
            packets: filteredPackets.slice(start, end),
            topSpacerHeight: start * PACKET_TRACE_ROW_HEIGHT,
            bottomSpacerHeight:
                Math.max(0, filteredPackets.length - end) *
                PACKET_TRACE_ROW_HEIGHT,
        };
    }, [filteredPackets, packetTraceScrollTop]);

    const getEffectiveOutputDir = async () => {
        if (outputDirectory) return outputDirectory;
        if (scenarioFilePath) {
            return dirname(scenarioFilePath);
        }
        return '.';
    };

    const handleStart = async () => {
        const scenarioState = useScenarioStore.getState();
        const validationMessage =
            getBlockingFmcwValidationMessage(scenarioState);
        if (validationMessage) {
            showError(`FMCW validation failed: ${validationMessage}`);
            return;
        }

        const errors = validateVita49Config(config);
        if (errors.length > 0) {
            showError(errors.join(' '));
            return;
        }

        const expected = deriveExpectedVita49Streams(scenarioState);
        if (expected.length === 0) {
            showWarning('No receiver streams are configured.');
        }

        setMetadataExportPath(null);
        startRun(expected);
        try {
            await invoke('set_output_directory', {
                dir: await getEffectiveOutputDir(),
            });
            await scenarioState.syncBackend();
            await invoke('start_vita49_stream', {
                config: toVita49BackendConfig(config),
            });
        } catch (err) {
            const message = err instanceof Error ? err.message : String(err);
            failRun(message);
            showError(`Failed to start VITA49 streaming: ${message}`);
        }
    };

    const handleStop = async () => {
        markStopping();
        try {
            await invoke('stop_simulation');
        } catch (err) {
            const message = err instanceof Error ? err.message : String(err);
            showError(`Failed to stop simulation: ${message}`);
        }
    };

    const exportMetadataJson = async () => {
        try {
            const outputPath = await invoke<string>(
                'export_output_metadata_json'
            );
            setMetadataExportPath(outputPath);
            showSuccess(`Metadata JSON saved to ${outputPath}`);
        } catch (err) {
            const message = err instanceof Error ? err.message : String(err);
            showError(`Failed to export metadata JSON: ${message}`);
        }
    };

    return (
        <Box sx={{ p: 3, height: '100%', overflowY: 'auto' }}>
            <Stack
                direction="row"
                spacing={2}
                alignItems="center"
                sx={{ mb: 2 }}
            >
                <Typography variant="h4">VITA49 Streams</Typography>
                <Chip
                    label={runState}
                    color={isRunning ? 'primary' : 'default'}
                />
            </Stack>

            {error && (
                <Alert severity="error" sx={{ mb: 2 }}>
                    {error}
                </Alert>
            )}

            <Grid container spacing={2} sx={{ mb: 2 }}>
                <Grid size={{ xs: 12, lg: 4 }}>
                    <Paper variant="outlined" sx={{ p: 2, height: '100%' }}>
                        <Typography variant="h6" sx={{ mb: 2 }}>
                            Runtime
                        </Typography>
                        <Stack spacing={2}>
                            <TextField
                                label="Host"
                                size="small"
                                value={config.host}
                                disabled={isRunning}
                                onChange={(event) =>
                                    setConfig({ host: event.target.value })
                                }
                            />
                            <TextField
                                label="Port"
                                size="small"
                                type="number"
                                value={config.port}
                                disabled={isRunning}
                                onChange={(event) =>
                                    setConfig({
                                        port: Number(event.target.value),
                                    })
                                }
                            />
                            <TextField
                                label="Full-scale"
                                size="small"
                                type="number"
                                value={config.fullscale}
                                disabled={isRunning}
                                onChange={(event) =>
                                    setConfig({
                                        fullscale: Number(event.target.value),
                                    })
                                }
                            />
                            <FormControl size="small">
                                <InputLabel id="vita49-epoch-mode-label">
                                    Epoch
                                </InputLabel>
                                <Select
                                    labelId="vita49-epoch-mode-label"
                                    label="Epoch"
                                    value={config.epochMode}
                                    disabled={isRunning}
                                    onChange={(event) =>
                                        setConfig({
                                            epochMode: event.target
                                                .value as typeof config.epochMode,
                                        })
                                    }
                                >
                                    <MenuItem value="auto">Auto</MenuItem>
                                    <MenuItem value="fixed">Fixed</MenuItem>
                                </Select>
                            </FormControl>
                            {config.epochMode === 'fixed' && (
                                <TextField
                                    label="Unix ns"
                                    size="small"
                                    value={config.epochUnixNanoseconds}
                                    disabled={isRunning}
                                    onChange={(event) =>
                                        setConfig({
                                            epochUnixNanoseconds:
                                                event.target.value,
                                        })
                                    }
                                />
                            )}
                            <TextField
                                label="Max UDP payload"
                                size="small"
                                type="number"
                                value={config.maxUdpPayload}
                                disabled={isRunning}
                                onChange={(event) =>
                                    setConfig({
                                        maxUdpPayload: Number(
                                            event.target.value
                                        ),
                                    })
                                }
                            />
                            <TextField
                                label="Queue depth"
                                size="small"
                                type="number"
                                value={config.queueDepth}
                                disabled={isRunning}
                                onChange={(event) =>
                                    setConfig({
                                        queueDepth: Number(event.target.value),
                                    })
                                }
                            />
                            <Stack
                                direction="row"
                                alignItems="center"
                                justifyContent="space-between"
                                spacing={1}
                            >
                                <Typography variant="body2">
                                    Packet trace
                                </Typography>
                                <Switch
                                    checked={config.traceEnabled}
                                    disabled={isRunning}
                                    onChange={(event) =>
                                        setConfig({
                                            traceEnabled: event.target.checked,
                                        })
                                    }
                                />
                            </Stack>
                            <TextField
                                label="Trace ring"
                                size="small"
                                type="number"
                                value={config.packetTraceRingSize}
                                disabled={isRunning}
                                onChange={(event) =>
                                    setConfig({
                                        packetTraceRingSize: Number(
                                            event.target.value
                                        ),
                                    })
                                }
                            />
                            {configErrors.length > 0 && (
                                <Alert severity="warning">
                                    {configErrors.join(' ')}
                                </Alert>
                            )}
                            <Stack direction="row" spacing={1}>
                                <Button
                                    variant="contained"
                                    startIcon={<PlayCircleOutlineIcon />}
                                    disabled={
                                        isRunning || configErrors.length > 0
                                    }
                                    onClick={handleStart}
                                >
                                    Start
                                </Button>
                                <Button
                                    variant="outlined"
                                    color="error"
                                    startIcon={<StopCircleIcon />}
                                    disabled={runState !== 'running'}
                                    onClick={handleStop}
                                >
                                    Stop
                                </Button>
                            </Stack>
                        </Stack>
                    </Paper>
                </Grid>

                <Grid size={{ xs: 12, lg: 8 }}>
                    <Stack spacing={2}>
                        <Paper variant="outlined" sx={{ p: 2 }}>
                            <Stack
                                direction={{ xs: 'column', md: 'row' }}
                                spacing={2}
                                justifyContent="space-between"
                            >
                                <Box>
                                    <Typography variant="overline">
                                        Endpoint
                                    </Typography>
                                    <Typography variant="h6">
                                        {config.host}:{config.port}
                                    </Typography>
                                </Box>
                                <Box>
                                    <Typography variant="overline">
                                        Profile
                                    </Typography>
                                    <Typography variant="h6">
                                        {finalVita49Metadata?.class_id ??
                                            '0xFA52530001000101'}
                                    </Typography>
                                </Box>
                                <Box>
                                    <Typography variant="overline">
                                        Epoch
                                    </Typography>
                                    <Typography variant="h6">
                                        {formatExact(
                                            streamStats?.epoch_unix_nanoseconds ??
                                                finalVita49Metadata?.epoch_unix_nanoseconds
                                        )}
                                    </Typography>
                                </Box>
                            </Stack>
                        </Paper>

                        <Grid container spacing={1}>
                            {[
                                { label: 'Packets', value: aggregate.packets },
                                { label: 'Samples', value: aggregate.samples },
                                { label: 'Drops', value: aggregate.drops },
                                {
                                    label: 'Late',
                                    value: totalLatePacketCount(
                                        aggregate.lateData,
                                        aggregate.lateContext
                                    ),
                                    lateInfo: true,
                                },
                                {
                                    label: 'Over-range',
                                    value: aggregate.overRange,
                                },
                                { label: 'Context', value: aggregate.context },
                            ].map(({ label, value, lateInfo }) => (
                                <Grid size={{ xs: 6, md: 2 }} key={label}>
                                    <Paper variant="outlined" sx={{ p: 1.5 }}>
                                        <Stack
                                            direction="row"
                                            alignItems="center"
                                        >
                                            <Typography variant="caption">
                                                {label}
                                            </Typography>
                                            {lateInfo && (
                                                <Vita49LatePacketInfo
                                                    dataPacketCount={
                                                        aggregate.lateData
                                                    }
                                                    contextPacketCount={
                                                        aggregate.lateContext
                                                    }
                                                />
                                            )}
                                        </Stack>
                                        <Typography variant="h6">
                                            {formatMetric(value)}
                                        </Typography>
                                    </Paper>
                                </Grid>
                            ))}
                        </Grid>
                    </Stack>
                </Grid>
            </Grid>

            <Paper variant="outlined" sx={{ p: 2, mb: 2 }}>
                <Stack
                    direction="row"
                    justifyContent="space-between"
                    alignItems="center"
                    sx={{ mb: 1 }}
                >
                    <Typography variant="h6">Streams</Typography>
                    <Button
                        variant="outlined"
                        startIcon={<SaveAltIcon />}
                        disabled={!finalVita49Metadata}
                        onClick={exportMetadataJson}
                    >
                        Export JSON
                    </Button>
                </Stack>
                {metadataExportPath && (
                    <Typography
                        variant="body2"
                        color="text.secondary"
                        sx={{ mb: 1, overflowWrap: 'anywhere' }}
                    >
                        {metadataExportPath}
                    </Typography>
                )}
                <TableContainer>
                    <Table size="small">
                        <TableHead>
                            <TableRow>
                                <TableCell>Receiver</TableCell>
                                <TableCell>Stream ID</TableCell>
                                <TableCell align="right">Rate</TableCell>
                                <TableCell align="right">RF</TableCell>
                                <TableCell align="right">Packets</TableCell>
                                <TableCell align="right">Samples</TableCell>
                                <TableCell align="right">Drops</TableCell>
                                <TableCell align="right">
                                    <Stack
                                        direction="row"
                                        alignItems="center"
                                        justifyContent="flex-end"
                                    >
                                        Late
                                        <Vita49LatePacketInfo
                                            dataPacketCount={aggregate.lateData}
                                            contextPacketCount={
                                                aggregate.lateContext
                                            }
                                        />
                                    </Stack>
                                </TableCell>
                                <TableCell align="right">Context</TableCell>
                                <TableCell>Simulation span</TableCell>
                                <TableCell>UTC span</TableCell>
                            </TableRow>
                        </TableHead>
                        <TableBody>
                            {streamRows.map((row) => (
                                <TableRow key={row.key}>
                                    <TableCell>
                                        <Stack spacing={0.25}>
                                            <Typography variant="body2">
                                                {row.receiverName}
                                            </Typography>
                                            <Typography
                                                variant="caption"
                                                color="text.secondary"
                                            >
                                                {row.platformName
                                                    ? `${row.platformName} / ${row.mode}`
                                                    : row.mode}
                                            </Typography>
                                        </Stack>
                                    </TableCell>
                                    <TableCell>
                                        {formatStreamId(row.streamId)}
                                    </TableCell>
                                    <TableCell align="right">
                                        {formatMetric(row.sampleRate)}
                                    </TableCell>
                                    <TableCell align="right">
                                        {formatMetric(row.referenceFrequency)}
                                    </TableCell>
                                    <TableCell align="right">
                                        {formatMetric(row.packetsEmitted)}
                                    </TableCell>
                                    <TableCell align="right">
                                        {formatMetric(row.samplesEmitted)}
                                    </TableCell>
                                    <TableCell align="right">
                                        {formatMetric(row.packetsDropped)}
                                    </TableCell>
                                    <TableCell align="right">
                                        <Stack
                                            direction="row"
                                            alignItems="center"
                                            justifyContent="flex-end"
                                        >
                                            {formatMetric(
                                                totalLatePacketCount(
                                                    row.lateDataPacketCount,
                                                    row.lateContextPacketCount
                                                )
                                            )}
                                            <Vita49LatePacketInfo
                                                dataPacketCount={
                                                    row.lateDataPacketCount
                                                }
                                                contextPacketCount={
                                                    row.lateContextPacketCount
                                                }
                                            />
                                        </Stack>
                                    </TableCell>
                                    <TableCell align="right">
                                        {formatMetric(row.contextPackets)}
                                    </TableCell>
                                    <TableCell>
                                        {formatSimulationSpan(
                                            row.firstSampleTime,
                                            row.endSampleTime
                                        )}
                                    </TableCell>
                                    <TableCell>
                                        {formatTimestampSpan(
                                            row.firstTimestamp,
                                            row.endTimestamp
                                        )}
                                    </TableCell>
                                </TableRow>
                            ))}
                            {streamRows.length === 0 && (
                                <TableRow>
                                    <TableCell colSpan={11}>
                                        <Typography color="text.secondary">
                                            No streams
                                        </Typography>
                                    </TableCell>
                                </TableRow>
                            )}
                        </TableBody>
                    </Table>
                </TableContainer>
            </Paper>

            <Paper variant="outlined" sx={{ p: 2 }}>
                <Stack
                    direction={{ xs: 'column', md: 'row' }}
                    spacing={2}
                    alignItems={{ xs: 'stretch', md: 'center' }}
                    sx={{ mb: 2 }}
                >
                    <Typography variant="h6" sx={{ flexGrow: 1 }}>
                        Packet Trace
                    </Typography>
                    <FilterListIcon color="action" />
                    <FormControl size="small" sx={{ minWidth: 160 }}>
                        <InputLabel id="vita49-stream-filter-label">
                            Stream
                        </InputLabel>
                        <Select
                            labelId="vita49-stream-filter-label"
                            label="Stream"
                            value={streamFilter}
                            onChange={(event) =>
                                setStreamFilter(event.target.value)
                            }
                        >
                            <MenuItem value="all">All</MenuItem>
                            {streamIdOptions.map((streamId) => (
                                <MenuItem
                                    value={String(streamId)}
                                    key={streamId}
                                >
                                    {formatStreamId(streamId)}
                                </MenuItem>
                            ))}
                        </Select>
                    </FormControl>
                    <FormControl size="small" sx={{ minWidth: 140 }}>
                        <InputLabel id="vita49-kind-filter-label">
                            Kind
                        </InputLabel>
                        <Select
                            labelId="vita49-kind-filter-label"
                            label="Kind"
                            value={packetKindFilter}
                            onChange={(event) =>
                                setPacketKindFilter(event.target.value)
                            }
                        >
                            <MenuItem value="all">All</MenuItem>
                            <MenuItem value="data">Data</MenuItem>
                            <MenuItem value="context">Context</MenuItem>
                        </Select>
                    </FormControl>
                    {[
                        ['Dropped', droppedOnly, setDroppedOnly],
                        ['Over-range', overRangeOnly, setOverRangeOnly],
                        ['Sample-loss', sampleLossOnly, setSampleLossOnly],
                    ].map(([label, checked, setter]) => (
                        <Stack
                            direction="row"
                            alignItems="center"
                            spacing={0.5}
                            key={label as string}
                        >
                            <Switch
                                size="small"
                                checked={checked as boolean}
                                onChange={(event) =>
                                    (setter as (value: boolean) => void)(
                                        event.target.checked
                                    )
                                }
                            />
                            <Typography variant="body2">
                                {label as string}
                            </Typography>
                        </Stack>
                    ))}
                </Stack>
                {omittedPacketTraceEvents > 0 && (
                    <Alert severity="info" sx={{ mb: 2 }}>
                        Showing last {formatMetric(packetTrace.length)} trace
                        events; {formatMetric(omittedPacketTraceEvents)} older
                        trace events discarded from trace history. Stream
                        packets and samples unaffected.
                    </Alert>
                )}
                <TableContainer
                    ref={packetTraceContainerRef}
                    onScroll={(event) =>
                        setPacketTraceScrollTop(event.currentTarget.scrollTop)
                    }
                    sx={{ maxHeight: PACKET_TRACE_TABLE_HEIGHT }}
                >
                    <Table size="small" stickyHeader>
                        <TableHead>
                            <TableRow>
                                <TableCell align="right">Seq</TableCell>
                                <TableCell>Event</TableCell>
                                <TableCell>Stream</TableCell>
                                <TableCell align="right">Bytes</TableCell>
                                <TableCell align="right">Samples</TableCell>
                                <TableCell align="right">t</TableCell>
                                <TableCell align="right">UTC</TableCell>
                                <TableCell>Flags</TableCell>
                            </TableRow>
                        </TableHead>
                        <TableBody>
                            {packetTraceWindow.topSpacerHeight > 0 && (
                                <TableRow
                                    sx={{
                                        height: `${packetTraceWindow.topSpacerHeight}px`,
                                    }}
                                >
                                    <TableCell
                                        colSpan={8}
                                        sx={{ p: 0, border: 0 }}
                                    />
                                </TableRow>
                            )}
                            {packetTraceWindow.packets.map((packet) => (
                                <TableRow key={packet.sequence}>
                                    <TableCell align="right">
                                        {packet.sequence}
                                    </TableCell>
                                    <TableCell>{packet.event}</TableCell>
                                    <TableCell>
                                        {formatStreamId(packet.stream_id)}
                                    </TableCell>
                                    <TableCell align="right">
                                        {formatMetric(packet.byte_count)}
                                    </TableCell>
                                    <TableCell align="right">
                                        {formatMetric(packet.sample_count)}
                                    </TableCell>
                                    <TableCell align="right">
                                        {formatSeconds(
                                            packet.first_sample_time
                                        )}
                                    </TableCell>
                                    <TableCell align="right">
                                        {formatVita49Timestamp(
                                            packet.timestamp
                                        )}
                                    </TableCell>
                                    <TableCell>
                                        <Stack direction="row" spacing={0.5}>
                                            {packet.dropped && (
                                                <Chip
                                                    label="drop"
                                                    size="small"
                                                    color="error"
                                                />
                                            )}
                                            {packet.over_range && (
                                                <Chip
                                                    label="over"
                                                    size="small"
                                                    color="warning"
                                                />
                                            )}
                                            {packet.sample_loss && (
                                                <Chip
                                                    label="loss"
                                                    size="small"
                                                    color="warning"
                                                />
                                            )}
                                        </Stack>
                                    </TableCell>
                                </TableRow>
                            ))}
                            {packetTraceWindow.bottomSpacerHeight > 0 && (
                                <TableRow
                                    sx={{
                                        height: `${packetTraceWindow.bottomSpacerHeight}px`,
                                    }}
                                >
                                    <TableCell
                                        colSpan={8}
                                        sx={{ p: 0, border: 0 }}
                                    />
                                </TableRow>
                            )}
                            {filteredPackets.length === 0 && (
                                <TableRow>
                                    <TableCell colSpan={8}>
                                        <Typography color="text.secondary">
                                            No packets
                                        </Typography>
                                    </TableCell>
                                </TableRow>
                            )}
                        </TableBody>
                    </Table>
                </TableContainer>
            </Paper>
        </Box>
    );
});
