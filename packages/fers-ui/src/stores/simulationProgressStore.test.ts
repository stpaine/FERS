// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import {
    normalizeSimulationOutputMetadata,
    type RawSimulationOutputMetadata,
    type SimulationOutputMetadata,
} from './simulationProgressStore';

const baseMetadata = {
    schema_version: 1,
    simulation_name: 'Metadata Fixture',
    output_directory: '/tmp/fers',
    start_time: 0,
    end_time: 1,
    sampling_rate: 1e6,
    oversample_ratio: 1,
} satisfies Omit<SimulationOutputMetadata, 'files'>;

describe('simulation output metadata', () => {
    test('keeps FMCW metadata and streaming segments', () => {
        const metadata: SimulationOutputMetadata = {
            ...baseMetadata,
            files: [
                {
                    receiver_id: 10,
                    receiver_name: 'FMCW Rx',
                    mode: 'fmcw',
                    path: '/tmp/fers/fmcw.h5',
                    total_samples: 128,
                    sample_start: 0,
                    sample_end_exclusive: 128,
                    pulse_count: 0,
                    min_pulse_length_samples: 0,
                    max_pulse_length_samples: 0,
                    uniform_pulse_length: true,
                    chunks: [],
                    streaming_segments: [
                        {
                            start_time: 0,
                            end_time: 0.001,
                            sample_count: 128,
                            sample_start: 0,
                            sample_end_exclusive: 128,
                            first_chirp_start_time: 0,
                            emitted_chirp_count: 4,
                        },
                    ],
                    fmcw: {
                        chirp_bandwidth: 20e6,
                        chirp_duration: 250e-6,
                        chirp_period: 500e-6,
                        chirp_rate: 80e9,
                        chirp_rate_signed: -80e9,
                        chirp_direction: 'down',
                        start_frequency_offset: 0,
                        chirp_count: 4,
                    },
                    fmcw_sources: [
                        {
                            transmitter_id: 11,
                            transmitter_name: 'FMCW Tx',
                            waveform_id: 12,
                            waveform_name: 'FMCW Wave',
                            carrier_frequency: 10e9,
                            chirp_bandwidth: 20e6,
                            chirp_duration: 250e-6,
                            chirp_period: 500e-6,
                            chirp_rate: 80e9,
                            chirp_rate_signed: -80e9,
                            chirp_direction: 'down',
                            start_frequency_offset: 0,
                            chirp_count: 4,
                            segments: [
                                {
                                    start_time: 0,
                                    end_time: 0.001,
                                    first_chirp_start_time: 0,
                                    emitted_chirp_count: 4,
                                },
                            ],
                        },
                    ],
                },
            ],
        };

        expect(
            normalizeSimulationOutputMetadata(metadata).files[0]
        ).toMatchObject({
            mode: 'fmcw',
            streaming_segments: [
                {
                    first_chirp_start_time: 0,
                    emitted_chirp_count: 4,
                },
            ],
            fmcw: {
                chirp_bandwidth: 20e6,
                chirp_duration: 250e-6,
                chirp_period: 500e-6,
                chirp_direction: 'down',
            },
            fmcw_sources: [
                {
                    transmitter_id: 11,
                    waveform_id: 12,
                    segments: [
                        {
                            first_chirp_start_time: 0,
                            emitted_chirp_count: 4,
                        },
                    ],
                },
            ],
        });
    });

    test('adapts legacy cw_segments metadata to streaming_segments', () => {
        const legacyMetadata = {
            ...baseMetadata,
            files: [
                {
                    receiver_id: 20,
                    receiver_name: 'CW Rx',
                    mode: 'cw',
                    path: '/tmp/fers/cw.h5',
                    total_samples: 64,
                    sample_start: 0,
                    sample_end_exclusive: 64,
                    pulse_count: 0,
                    min_pulse_length_samples: 0,
                    max_pulse_length_samples: 0,
                    uniform_pulse_length: true,
                    chunks: [],
                    cw_segments: [
                        {
                            start_time: 0,
                            end_time: 0.001,
                            sample_count: 64,
                            sample_start: 0,
                            sample_end_exclusive: 64,
                        },
                    ],
                },
            ],
        } satisfies RawSimulationOutputMetadata;

        const normalizedFile =
            normalizeSimulationOutputMetadata(legacyMetadata).files[0];

        expect(normalizedFile.streaming_segments).toHaveLength(1);
        expect(normalizedFile.fmcw_sources).toEqual([]);
        expect('cw_segments' in normalizedFile).toBe(false);
    });
});
