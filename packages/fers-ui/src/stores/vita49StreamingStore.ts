// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { create } from 'zustand';
import { persist } from 'zustand/middleware';
import type { ScenarioData } from './scenarioStore';
import type {
    SimulationOutputMetadata,
    SimulationOutputVita49Metadata,
    SimulationOutputVita49Timestamp,
} from './simulationProgressStore';

export type Vita49Timestamp = SimulationOutputVita49Timestamp;

export type Vita49RunState =
    | 'idle'
    | 'running'
    | 'stopping'
    | 'draining'
    | 'completed'
    | 'failed'
    | 'cancelled';

export type Vita49EpochMode = 'auto' | 'fixed';

export type Vita49RuntimeConfig = {
    host: string;
    port: number;
    fullscale: number;
    epochMode: Vita49EpochMode;
    epochUnixNanoseconds: string;
    maxUdpPayload: number;
    queueDepth: number;
    traceEnabled: boolean;
    packetTraceRingSize: number;
};

export type Vita49BackendConfig = {
    host: string;
    port: number;
    fullscale: number;
    epoch_unix_nanoseconds: string | null;
    max_udp_payload: number;
    queue_depth: number;
    trace_enabled: boolean;
    packet_trace_ring_size: number;
};

export type Vita49StreamCounter = {
    receiver_id: number;
    receiver_name: string;
    stream_id: number;
    mode?: string | null;
    sample_rate: number;
    reference_frequency: number;
    packets_emitted: number;
    context_packets?: number;
    context_packet_count?: number;
    samples_emitted: number;
    packets_dropped: number;
    samples_dropped: number;
    over_range_count: number;
    late_packet_count: number;
    first_sample_time: number | null;
    end_sample_time: number | null;
    first_timestamp: Vita49Timestamp | null;
    end_timestamp: Vita49Timestamp | null;
};

export type Vita49StreamStatsEvent = {
    mode: 'vita49_udp' | 'hdf5';
    epoch_unix_nanoseconds: string | null;
    streams: Vita49StreamCounter[];
};

export type Vita49PacketTraceEvent = {
    sequence: number;
    event: 'data' | 'context' | 'drop' | 'failure' | string;
    stream_id: number;
    byte_count: number;
    sample_count: number;
    first_sample_time: number;
    timestamp: Vita49Timestamp | null;
    data_packet: boolean;
    context_packet: boolean;
    dropped: boolean;
    over_range: boolean;
    sample_loss: boolean;
};

export type Vita49TelemetryPoll = {
    stats: Vita49StreamStatsEvent | null;
    packets: Vita49PacketTraceEvent[];
    omitted_packet_trace_events: number;
    has_more: boolean;
};

export type Vita49StreamRow = {
    key: string;
    receiverId: number;
    receiverName: string;
    platformName: string;
    mode: Vita49StreamMode;
    streamId: number | null;
    sampleRate: number | null;
    referenceFrequency: number | null;
    packetsEmitted: number;
    contextPackets: number;
    samplesEmitted: number;
    packetsDropped: number;
    samplesDropped: number;
    overRangeCount: number;
    latePacketCount: number;
    firstSampleTime: number | null;
    endSampleTime: number | null;
    firstTimestamp: Vita49Timestamp | null;
    endTimestamp: Vita49Timestamp | null;
    backendObserved: boolean;
};

export type Vita49StreamMode = 'pulsed' | 'cw' | 'fmcw' | 'sfcw' | 'unknown';

export const DEFAULT_VITA49_CONFIG: Vita49RuntimeConfig = {
    host: '127.0.0.1',
    port: 4991,
    fullscale: 1.0,
    epochMode: 'auto',
    epochUnixNanoseconds: '',
    maxUdpPayload: 1400,
    queueDepth: 1024,
    traceEnabled: true,
    packetTraceRingSize: 500,
};

const MAX_VRT_EPOCH_NS = 4_294_967_295_999_999_999n;
const MAX_VITA49_QUEUE_DEPTH = 4_294_967_295;

const isPositiveFinite = (value: number) => Number.isFinite(value) && value > 0;

const isIntegerInRange = (value: number, min: number, max: number) =>
    Number.isInteger(value) && value >= min && value <= max;

export const validateVita49Config = (config: Vita49RuntimeConfig): string[] => {
    const errors: string[] = [];

    if (config.host.trim().length === 0) {
        errors.push('Host is required.');
    }
    if (!isIntegerInRange(config.port, 1, 65535)) {
        errors.push('Port must be an integer from 1 to 65535.');
    }
    if (!isPositiveFinite(config.fullscale)) {
        errors.push('Full-scale must be a positive finite number.');
    }
    if (!isIntegerInRange(config.maxUdpPayload, 64, 65507)) {
        errors.push('Max UDP payload must be an integer from 64 to 65507.');
    }
    if (!isIntegerInRange(config.queueDepth, 1, MAX_VITA49_QUEUE_DEPTH)) {
        errors.push('Queue depth must be an integer from 1 to 4294967295.');
    }
    if (
        !isIntegerInRange(
            config.packetTraceRingSize,
            1,
            Number.MAX_SAFE_INTEGER
        )
    ) {
        errors.push('Packet trace ring size must be a positive integer.');
    }
    if (config.epochMode === 'fixed') {
        try {
            if (!/^\d+$/.test(config.epochUnixNanoseconds.trim())) {
                throw new Error('not an integer');
            }
            const epoch = BigInt(config.epochUnixNanoseconds.trim());
            if (epoch > MAX_VRT_EPOCH_NS) {
                errors.push(
                    'Fixed epoch must fit the VRT 32-bit UTC seconds field.'
                );
            }
        } catch {
            errors.push('Fixed epoch must be a Unix nanosecond integer.');
        }
    }

    return errors;
};

export const toVita49BackendConfig = (
    config: Vita49RuntimeConfig
): Vita49BackendConfig => {
    const errors = validateVita49Config(config);
    if (errors.length > 0) {
        throw new Error(errors.join(' '));
    }

    return {
        host: config.host.trim(),
        port: config.port,
        fullscale: config.fullscale,
        epoch_unix_nanoseconds:
            config.epochMode === 'fixed'
                ? config.epochUnixNanoseconds.trim()
                : null,
        max_udp_payload: config.maxUdpPayload,
        queue_depth: config.queueDepth,
        trace_enabled: config.traceEnabled,
        packet_trace_ring_size: config.packetTraceRingSize,
    };
};

const normalizeStreamMode = (
    mode: string | null | undefined
): Vita49StreamMode =>
    mode === 'pulsed' || mode === 'cw' || mode === 'fmcw' || mode === 'sfcw'
        ? mode
        : 'unknown';

export const deriveExpectedVita49Streams = (
    scenario: Pick<ScenarioData, 'globalParameters' | 'platforms' | 'waveforms'>
): Vita49StreamRow[] => {
    const waveformById = new Map(
        scenario.waveforms.map((waveform) => [waveform.id, waveform])
    );
    const rows: Vita49StreamRow[] = [];

    for (const platform of scenario.platforms) {
        for (const component of platform.components) {
            if (
                component.type !== 'receiver' &&
                component.type !== 'monostatic'
            ) {
                continue;
            }
            const waveform =
                component.type === 'monostatic' && component.waveformId
                    ? waveformById.get(component.waveformId)
                    : undefined;
            const receiverId =
                component.type === 'monostatic'
                    ? Number(component.rxId)
                    : Number(component.id);

            rows.push({
                key: `${receiverId}:${component.radarType}`,
                receiverId,
                receiverName: component.name,
                platformName: platform.name,
                mode: normalizeStreamMode(component.radarType),
                streamId: null,
                sampleRate:
                    component.fmcwModeConfig?.if_sample_rate ??
                    scenario.globalParameters.rate,
                referenceFrequency: waveform?.carrier_frequency ?? null,
                packetsEmitted: 0,
                contextPackets: 0,
                samplesEmitted: 0,
                packetsDropped: 0,
                samplesDropped: 0,
                overRangeCount: 0,
                latePacketCount: 0,
                firstSampleTime: null,
                endSampleTime: null,
                firstTimestamp: null,
                endTimestamp: null,
                backendObserved: false,
            });
        }
    }

    return rows;
};

const rowFromCounter = (counter: Vita49StreamCounter): Vita49StreamRow => ({
    key: `${counter.receiver_id}:${normalizeStreamMode(counter.mode)}:${counter.stream_id}`,
    receiverId: counter.receiver_id,
    receiverName: counter.receiver_name,
    platformName: '',
    mode: normalizeStreamMode(counter.mode),
    streamId: counter.stream_id,
    sampleRate: counter.sample_rate,
    referenceFrequency: counter.reference_frequency,
    packetsEmitted: counter.packets_emitted,
    contextPackets:
        counter.context_packets ?? counter.context_packet_count ?? 0,
    samplesEmitted: counter.samples_emitted,
    packetsDropped: counter.packets_dropped,
    samplesDropped: counter.samples_dropped,
    overRangeCount: counter.over_range_count,
    latePacketCount: counter.late_packet_count,
    firstSampleTime: counter.first_sample_time,
    endSampleTime: counter.end_sample_time,
    firstTimestamp: counter.first_timestamp,
    endTimestamp: counter.end_timestamp,
    backendObserved: true,
});

export const mergeVita49StreamRows = (
    expectedRows: Vita49StreamRow[],
    stats: Vita49StreamStatsEvent | null
): Vita49StreamRow[] => {
    if (!stats) {
        return expectedRows;
    }

    const rowsByKey = new Map(expectedRows.map((row) => [row.key, { ...row }]));
    const expectedByReceiver = new Map<number, Vita49StreamRow[]>();
    for (const row of expectedRows) {
        expectedByReceiver.set(row.receiverId, [
            ...(expectedByReceiver.get(row.receiverId) ?? []),
            row,
        ]);
    }

    const findExpectedRow = (counter: Vita49StreamCounter) => {
        const observedMode = normalizeStreamMode(counter.mode);
        if (observedMode !== 'unknown') {
            const exact = rowsByKey.get(
                `${counter.receiver_id}:${observedMode}`
            );
            if (exact) return exact;
        }
        const receiverRows = expectedByReceiver.get(counter.receiver_id) ?? [];
        return receiverRows.length === 1
            ? rowsByKey.get(receiverRows[0].key)
            : undefined;
    };

    const rowsByResolvedKey = new Map(
        expectedRows.map((row) => [row.key, { ...row }])
    );

    for (const counter of stats.streams) {
        const existing = findExpectedRow(counter);
        const observed = rowFromCounter(counter);
        const modeFromBackend = normalizeStreamMode(counter.mode);
        const key = existing?.key ?? observed.key;
        rowsByResolvedKey.set(key, {
            ...(existing ?? observed),
            receiverName: counter.receiver_name || existing?.receiverName || '',
            mode:
                modeFromBackend !== 'unknown'
                    ? modeFromBackend
                    : (existing?.mode ?? observed.mode),
            streamId: counter.stream_id,
            sampleRate: counter.sample_rate,
            referenceFrequency: counter.reference_frequency,
            packetsEmitted: counter.packets_emitted,
            contextPackets:
                counter.context_packets ?? counter.context_packet_count ?? 0,
            samplesEmitted: counter.samples_emitted,
            packetsDropped: counter.packets_dropped,
            samplesDropped: counter.samples_dropped,
            overRangeCount: counter.over_range_count,
            latePacketCount: counter.late_packet_count,
            firstSampleTime: counter.first_sample_time,
            endSampleTime: counter.end_sample_time,
            firstTimestamp: counter.first_timestamp,
            endTimestamp: counter.end_timestamp,
            backendObserved: true,
        });
    }

    const modeOrder: Record<Vita49StreamMode, number> = {
        pulsed: 0,
        cw: 1,
        fmcw: 2,
        sfcw: 3,
        unknown: 4,
    };
    return Array.from(rowsByResolvedKey.values()).sort(
        (a, b) =>
            a.receiverId - b.receiverId ||
            modeOrder[a.mode] - modeOrder[b.mode] ||
            (a.streamId ?? 0) - (b.streamId ?? 0)
    );
};

type Vita49StreamingStore = {
    config: Vita49RuntimeConfig;
    runState: Vita49RunState;
    expectedStreams: Vita49StreamRow[];
    streamStats: Vita49StreamStatsEvent | null;
    packetTrace: Vita49PacketTraceEvent[];
    omittedPacketTraceEvents: number;
    finalMetadata: SimulationOutputMetadata | null;
    finalVita49Metadata: SimulationOutputVita49Metadata | null;
    error: string | null;
    setConfig: (config: Partial<Vita49RuntimeConfig>) => void;
    startRun: (expectedStreams: Vita49StreamRow[]) => void;
    markStopping: () => void;
    markDraining: () => void;
    setStreamStats: (stats: Vita49StreamStatsEvent) => void;
    appendPacketBatch: (
        packets: Vita49PacketTraceEvent[],
        omittedPacketTraceEvents?: number
    ) => void;
    completeRun: (metadata: SimulationOutputMetadata | null) => void;
    cancelRun: (metadata: SimulationOutputMetadata | null) => void;
    failRun: (error: string) => void;
    resetTrace: () => void;
};

export const useVita49StreamingStore = create<Vita49StreamingStore>()(
    persist(
        (set, get) => ({
            config: DEFAULT_VITA49_CONFIG,
            runState: 'idle',
            expectedStreams: [],
            streamStats: null,
            packetTrace: [],
            omittedPacketTraceEvents: 0,
            finalMetadata: null,
            finalVita49Metadata: null,
            error: null,

            setConfig: (config) =>
                set((state) => ({ config: { ...state.config, ...config } })),
            startRun: (expectedStreams) =>
                set({
                    runState: 'running',
                    expectedStreams,
                    streamStats: null,
                    packetTrace: [],
                    omittedPacketTraceEvents: 0,
                    finalMetadata: null,
                    finalVita49Metadata: null,
                    error: null,
                }),
            markStopping: () =>
                set((state) => ({
                    runState:
                        state.runState === 'running'
                            ? 'stopping'
                            : state.runState,
                })),
            markDraining: () =>
                set((state) => ({
                    runState:
                        state.runState === 'running' ||
                        state.runState === 'stopping'
                            ? 'draining'
                            : state.runState,
                })),
            setStreamStats: (streamStats) => set({ streamStats }),
            appendPacketBatch: (packets, omittedPacketTraceEvents = 0) =>
                set((state) => {
                    const ringSize = state.config.packetTraceRingSize;
                    const combined = [...state.packetTrace, ...packets];
                    const overflow = Math.max(0, combined.length - ringSize);
                    return {
                        packetTrace:
                            overflow > 0 ? combined.slice(overflow) : combined,
                        omittedPacketTraceEvents:
                            state.omittedPacketTraceEvents +
                            omittedPacketTraceEvents +
                            overflow,
                    };
                }),
            completeRun: (metadata) =>
                set({
                    runState: 'completed',
                    finalMetadata: metadata,
                    finalVita49Metadata: metadata?.vita49 ?? null,
                    error: null,
                }),
            cancelRun: (metadata) =>
                set({
                    runState: 'cancelled',
                    finalMetadata: metadata,
                    finalVita49Metadata: metadata?.vita49 ?? null,
                    error: null,
                }),
            failRun: (error) =>
                set({
                    runState: 'failed',
                    error,
                }),
            resetTrace: () =>
                set({ packetTrace: [], omittedPacketTraceEvents: 0 }),
        }),
        {
            name: 'fers-vita49-streaming',
            partialize: (state) => ({ config: state.config }),
            merge: (persisted, current) => {
                const persistedState =
                    persisted && typeof persisted === 'object'
                        ? (persisted as Partial<
                              Pick<Vita49StreamingStore, 'config'>
                          >)
                        : {};
                return {
                    ...current,
                    ...persistedState,
                    config: {
                        ...DEFAULT_VITA49_CONFIG,
                        ...persistedState.config,
                    },
                };
            },
        }
    )
);
