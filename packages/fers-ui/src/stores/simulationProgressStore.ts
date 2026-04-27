// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { create } from 'zustand';

export type SimulationProgressState = {
    message: string;
    current: number;
    total: number;
    details?: SimulationProgressDetail[];
};

export type SimulationProgressDetail = {
    id: string;
    message: string;
    current: number;
    total: number;
};

export type SimulationRunStatus = 'idle' | 'running' | 'completed' | 'failed';

export type SimulationOutputMode = 'pulsed' | 'cw' | 'fmcw';

export type SimulationOutputChunkMetadata = {
    chunk_index: number;
    i_dataset: string;
    q_dataset: string;
    start_time: number;
    sample_count: number;
    sample_start: number;
    sample_end_exclusive: number;
};

export type SimulationOutputStreamingSegmentMetadata = {
    start_time: number;
    end_time: number;
    sample_count: number;
    sample_start: number;
    sample_end_exclusive: number;
    first_chirp_start_time?: number;
    emitted_chirp_count?: number;
};

export type SimulationOutputFmcwMetadata = {
    chirp_bandwidth: number;
    chirp_duration: number;
    chirp_period: number;
    chirp_rate: number;
    start_frequency_offset: number;
    chirp_count?: number;
};

export type SimulationOutputFileMetadata = {
    receiver_id: number;
    receiver_name: string;
    mode: SimulationOutputMode;
    path: string;
    total_samples: number;
    sample_start: number;
    sample_end_exclusive: number;
    pulse_count: number;
    min_pulse_length_samples: number;
    max_pulse_length_samples: number;
    uniform_pulse_length: boolean;
    chunks: SimulationOutputChunkMetadata[];
    streaming_segments: SimulationOutputStreamingSegmentMetadata[];
    fmcw?: SimulationOutputFmcwMetadata;
};

export type SimulationOutputMetadata = {
    schema_version: number;
    simulation_name: string;
    output_directory: string;
    start_time: number;
    end_time: number;
    sampling_rate: number;
    oversample_ratio: number;
    files: SimulationOutputFileMetadata[];
};

export type RawSimulationOutputFileMetadata = Omit<
    SimulationOutputFileMetadata,
    'streaming_segments'
> & {
    streaming_segments?: SimulationOutputStreamingSegmentMetadata[];
    cw_segments?: SimulationOutputStreamingSegmentMetadata[];
};

export type RawSimulationOutputMetadata = Omit<
    SimulationOutputMetadata,
    'files'
> & {
    files: RawSimulationOutputFileMetadata[];
};

export const normalizeSimulationOutputMetadata = (
    metadata: RawSimulationOutputMetadata
): SimulationOutputMetadata => ({
    ...metadata,
    files: metadata.files.map((file) => {
        const { cw_segments, streaming_segments, ...normalizedFile } = file;
        return {
            ...normalizedFile,
            streaming_segments: streaming_segments ?? cw_segments ?? [],
        };
    }),
});

type SimulationProgressStore = {
    isSimulating: boolean;
    isGeneratingKml: boolean;
    simulationProgress: Record<string, SimulationProgressState>;
    simulationRunStatus: SimulationRunStatus;
    simulationRunError: string | null;
    simulationOutputMetadata: SimulationOutputMetadata | null;

    setIsSimulating: (isSimulating: boolean) => void;
    setIsGeneratingKml: (isGeneratingKml: boolean) => void;
    startSimulationRun: () => void;
    setSimulationProgressSnapshot: (
        progress: Record<string, SimulationProgressState>
    ) => void;
    setSimulationOutputMetadata: (
        metadata: SimulationOutputMetadata | null
    ) => void;
    completeSimulationRun: () => void;
    failSimulationRun: (errorMessage: string) => void;
    clearSimulationProgress: () => void;
};

export const useSimulationProgressStore = create<SimulationProgressStore>()(
    (set) => ({
        isSimulating: false,
        isGeneratingKml: false,
        simulationProgress: {},
        simulationRunStatus: 'idle',
        simulationRunError: null,
        simulationOutputMetadata: null,

        setIsSimulating: (isSimulating) => set({ isSimulating }),
        setIsGeneratingKml: (isGeneratingKml) => set({ isGeneratingKml }),
        startSimulationRun: () =>
            set({
                isSimulating: true,
                simulationProgress: {},
                simulationRunStatus: 'running',
                simulationRunError: null,
                simulationOutputMetadata: null,
            }),
        setSimulationProgressSnapshot: (progress) =>
            set({ simulationProgress: progress }),
        setSimulationOutputMetadata: (metadata) =>
            set({ simulationOutputMetadata: metadata }),
        completeSimulationRun: () =>
            set({
                isSimulating: false,
                simulationRunStatus: 'completed',
                simulationRunError: null,
            }),
        failSimulationRun: (errorMessage) =>
            set({
                isSimulating: false,
                simulationRunStatus: 'failed',
                simulationRunError: errorMessage,
            }),
        clearSimulationProgress: () =>
            set({
                simulationProgress: {},
                simulationRunStatus: 'idle',
                simulationRunError: null,
                simulationOutputMetadata: null,
            }),
    })
);
