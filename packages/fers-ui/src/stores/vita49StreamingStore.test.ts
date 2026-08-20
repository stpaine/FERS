// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import { defaultGlobalParameters } from './scenarioStore/defaults';
import type { ScenarioData } from './scenarioStore/types';
import {
    DEFAULT_VITA49_CONFIG,
    deriveExpectedVita49Streams,
    mergeVita49StreamRows,
    toVita49BackendConfig,
    useVita49StreamingStore,
    type Vita49PacketTraceEvent,
    validateVita49Config,
} from './vita49StreamingStore';

const packet = (sequence: number): Vita49PacketTraceEvent => ({
    sequence,
    event: 'data',
    stream_id: 10,
    byte_count: 128,
    sample_count: 24,
    first_sample_time: 0,
    timestamp: { integer_seconds: 1_700_000_000, fractional_picoseconds: 0 },
    data_packet: true,
    context_packet: false,
    dropped: false,
    over_range: false,
    sample_loss: false,
});

describe('VITA49 streaming config', () => {
    test('validates malformed endpoint and runtime bounds', () => {
        expect(
            validateVita49Config({
                ...DEFAULT_VITA49_CONFIG,
                host: ' ',
                port: 0,
                fullscale: 0,
                maxUdpPayload: 63,
                queueDepth: 0,
                packetTraceRingSize: 0,
            })
        ).toEqual([
            'Host is required.',
            'Port must be an integer from 1 to 65535.',
            'Full-scale must be a positive finite number.',
            'Max UDP payload must be an integer from 64 to 65507.',
            'Queue depth must be an integer from 1 to 4294967295.',
            'Packet trace ring size must be a positive integer.',
        ]);
    });

    test('rejects queue depth above uint32 max', () => {
        expect(
            validateVita49Config({
                ...DEFAULT_VITA49_CONFIG,
                queueDepth: 4_294_967_296,
            })
        ).toContain('Queue depth must be an integer from 1 to 4294967295.');
        expect(
            validateVita49Config({
                ...DEFAULT_VITA49_CONFIG,
                queueDepth: 4_294_967_295,
            })
        ).not.toContain('Queue depth must be an integer from 1 to 4294967295.');
    });

    test('validates fixed epoch and preserves ns as a string for Rust', () => {
        expect(
            validateVita49Config({
                ...DEFAULT_VITA49_CONFIG,
                epochMode: 'fixed',
                epochUnixNanoseconds: '4294967296000000000',
            })
        ).toContain('Fixed epoch must fit the VRT 32-bit UTC seconds field.');

        expect(
            toVita49BackendConfig({
                ...DEFAULT_VITA49_CONFIG,
                epochMode: 'fixed',
                epochUnixNanoseconds: '1700000000123456789',
            })
        ).toMatchObject({
            epoch_unix_nanoseconds: '1700000000123456789',
            trace_enabled: true,
            packet_trace_ring_size: DEFAULT_VITA49_CONFIG.packetTraceRingSize,
        });
    });
});

describe('VITA49 expected stream derivation', () => {
    test('covers pulsed, CW, and FMCW receiver streams', () => {
        const scenario: Pick<
            ScenarioData,
            'globalParameters' | 'platforms' | 'waveforms'
        > = {
            globalParameters: { ...defaultGlobalParameters, rate: 2e6 },
            waveforms: [
                {
                    id: '100',
                    type: 'Waveform',
                    name: 'CW',
                    waveformType: 'cw',
                    power: 1,
                    carrier_frequency: 9.6e9,
                },
            ],
            platforms: [
                {
                    id: '1',
                    type: 'Platform',
                    name: 'Platform A',
                    motionPath: {
                        interpolation: 'static',
                        waypoints: [
                            { id: '1', x: 0, y: 0, altitude: 0, time: 0 },
                        ],
                    },
                    rotation: {
                        type: 'fixed',
                        startAzimuth: 0,
                        startElevation: 0,
                        azimuthRate: 0,
                        elevationRate: 0,
                    },
                    components: [
                        {
                            id: '10',
                            type: 'receiver',
                            name: 'Rx',
                            radarType: 'pulsed',
                            window_skip: 1e-6,
                            window_length: 10e-6,
                            prf: 1000,
                            antennaId: null,
                            timingId: null,
                            noiseTemperature: null,
                            noDirectPaths: false,
                            noPropagationLoss: false,
                            schedule: [],
                        },
                        {
                            id: '11',
                            type: 'receiver',
                            name: 'CwRx',
                            radarType: 'cw',
                            window_skip: null,
                            window_length: null,
                            prf: null,
                            antennaId: null,
                            timingId: null,
                            noiseTemperature: null,
                            noDirectPaths: false,
                            noPropagationLoss: false,
                            schedule: [],
                        },
                        {
                            id: '12',
                            type: 'monostatic',
                            name: 'Mono',
                            txId: '13',
                            rxId: '14',
                            radarType: 'fmcw',
                            window_skip: null,
                            window_length: null,
                            prf: null,
                            antennaId: null,
                            waveformId: '100',
                            timingId: null,
                            noiseTemperature: null,
                            noDirectPaths: false,
                            noPropagationLoss: false,
                            fmcwModeConfig: { if_sample_rate: 250000 },
                            schedule: [],
                        },
                    ],
                },
            ],
        };

        expect(deriveExpectedVita49Streams(scenario)).toMatchObject([
            {
                receiverId: 10,
                receiverName: 'Rx',
                mode: 'pulsed',
                sampleRate: 2e6,
                referenceFrequency: null,
            },
            {
                receiverId: 11,
                receiverName: 'CwRx',
                mode: 'cw',
                sampleRate: 2e6,
                referenceFrequency: null,
            },
            {
                receiverId: 14,
                receiverName: 'Mono',
                mode: 'fmcw',
                sampleRate: 250000,
                referenceFrequency: 9.6e9,
            },
        ]);
    });
});

describe('VITA49 run lifecycle', () => {
    test('tracks draining between active streaming and completion', () => {
        useVita49StreamingStore.setState({
            runState: 'idle',
            expectedStreams: [],
            streamStats: null,
            packetTrace: [],
            omittedPacketTraceEvents: 0,
            finalMetadata: null,
            finalVita49Metadata: null,
            error: null,
        });

        useVita49StreamingStore.getState().startRun([]);
        expect(useVita49StreamingStore.getState().runState).toBe('running');

        useVita49StreamingStore.getState().markDraining();
        expect(useVita49StreamingStore.getState().runState).toBe('draining');

        useVita49StreamingStore.getState().completeRun(null);
        expect(useVita49StreamingStore.getState().runState).toBe('completed');
    });
});

describe('VITA49 telemetry rows', () => {
    test('maps exact VRT timestamp objects and simulation span', () => {
        const rows = mergeVita49StreamRows([], {
            mode: 'vita49_udp',
            epoch_unix_nanoseconds: '1700000000123456789',
            streams: [
                {
                    receiver_id: 7,
                    receiver_name: 'Rx',
                    stream_id: 1234,
                    mode: 'pulsed',
                    sample_rate: 100000,
                    reference_frequency: 10000000,
                    packets_emitted: 2930,
                    samples_emitted: 1000000,
                    packets_dropped: 0,
                    samples_dropped: 0,
                    over_range_count: 0,
                    late_data_packet_count: 3,
                    late_context_packet_count: 9,
                    context_packets: 12,
                    first_sample_time: 0,
                    end_sample_time: 10,
                    first_timestamp: {
                        integer_seconds: 1700000000,
                        fractional_picoseconds: 123456789000,
                    },
                    end_timestamp: {
                        integer_seconds: 1700000010,
                        fractional_picoseconds: 123456789000,
                    },
                },
            ],
        });

        expect(rows).toMatchObject([
            {
                receiverId: 7,
                mode: 'pulsed',
                samplesEmitted: 1000000,
                lateDataPacketCount: 3,
                lateContextPacketCount: 9,
                firstSampleTime: 0,
                endSampleTime: 10,
                firstTimestamp: {
                    integer_seconds: 1700000000,
                    fractional_picoseconds: 123456789000,
                },
                endTimestamp: {
                    integer_seconds: 1700000010,
                    fractional_picoseconds: 123456789000,
                },
            },
        ]);
    });

    test('does not default backend-only streams to CW', () => {
        const rows = mergeVita49StreamRows([], {
            mode: 'vita49_udp',
            epoch_unix_nanoseconds: null,
            streams: [
                {
                    receiver_id: 8,
                    receiver_name: 'BackendRx',
                    stream_id: 5678,
                    sample_rate: 100000,
                    reference_frequency: 10000000,
                    packets_emitted: 1,
                    samples_emitted: 10,
                    packets_dropped: 0,
                    samples_dropped: 0,
                    over_range_count: 0,
                    late_data_packet_count: 0,
                    late_context_packet_count: 0,
                    context_packets: 2,
                    first_sample_time: null,
                    end_sample_time: null,
                    first_timestamp: null,
                    end_timestamp: null,
                },
            ],
        });

        expect(rows).toMatchObject([
            {
                receiverId: 8,
                mode: 'unknown',
                backendObserved: true,
            },
        ]);
    });
});

describe('VITA49 packet trace ring', () => {
    test('bounds packets and counts omitted events', () => {
        useVita49StreamingStore.setState({
            config: { ...DEFAULT_VITA49_CONFIG, packetTraceRingSize: 3 },
            packetTrace: [],
            omittedPacketTraceEvents: 0,
        });

        useVita49StreamingStore
            .getState()
            .appendPacketBatch([1, 2, 3, 4, 5].map(packet));

        const state = useVita49StreamingStore.getState();
        expect(state.packetTrace.map((entry) => entry.sequence)).toEqual([
            3, 4, 5,
        ]);
        expect(state.omittedPacketTraceEvents).toBe(2);
    });

    test('adds backend omitted packet count when appending polled batches', () => {
        useVita49StreamingStore.setState({
            config: { ...DEFAULT_VITA49_CONFIG, packetTraceRingSize: 5 },
            packetTrace: [],
            omittedPacketTraceEvents: 0,
        });

        useVita49StreamingStore
            .getState()
            .appendPacketBatch([1, 2].map(packet), 7);

        const state = useVita49StreamingStore.getState();
        expect(state.packetTrace.map((entry) => entry.sequence)).toEqual([
            1, 2,
        ]);
        expect(state.omittedPacketTraceEvents).toBe(7);
    });
});
