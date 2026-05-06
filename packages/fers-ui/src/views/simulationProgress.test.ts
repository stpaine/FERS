// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import type { SimulationProgressState } from '@/stores/simulationProgressStore';
import {
    addSimulationProgressEvent,
    getSimulationProgressPercent,
    normalizeCompletedProgressSnapshot,
} from './simulationProgress';

const trackProgress = (events: SimulationProgressState[]) => {
    const tracked: Record<string, SimulationProgressState> = {};

    for (const event of events) {
        Object.assign(tracked, addSimulationProgressEvent(tracked, event));
    }

    return tracked;
};

describe('simulation progress normalization', () => {
    test('keeps CW receiver phases on one stable row', () => {
        const tracked = trackProgress([
            {
                message: 'Finalizing CW Receiver CWRadar',
                current: 0,
                total: 100,
            },
            {
                message: 'Rendering Interference for CWRadar',
                current: 25,
                total: 100,
            },
            {
                message: 'Applying Noise for CWRadar',
                current: 50,
                total: 100,
            },
            {
                message: 'Writing HDF5 for CWRadar',
                current: 75,
                total: 100,
            },
            {
                message: 'Finalized CWRadar',
                current: 100,
                total: 100,
            },
        ]);

        expect(Object.keys(tracked)).toEqual(['receiver:CWRadar']);
        expect(tracked['receiver:CWRadar']).toEqual({
            message: 'Finalized CWRadar',
            current: 100,
            total: 100,
            details: [
                {
                    id: 'streaming-finalizing',
                    message: 'Finalizing CW Receiver CWRadar',
                    current: 0,
                    total: 100,
                },
                {
                    id: 'streaming-rendering-interference',
                    message: 'Rendering Interference for CWRadar',
                    current: 25,
                    total: 100,
                },
                {
                    id: 'streaming-applying-noise',
                    message: 'Applying Noise for CWRadar',
                    current: 50,
                    total: 100,
                },
                {
                    id: 'streaming-writing-hdf5',
                    message: 'Writing HDF5 for CWRadar',
                    current: 75,
                    total: 100,
                },
                {
                    id: 'streaming-finalized',
                    message: 'Finalized CWRadar',
                    current: 100,
                    total: 100,
                },
            ],
        });
    });

    test('keeps FMCW receiver phases on one stable row', () => {
        const tracked = trackProgress([
            {
                message: 'Finalizing FMCW Receiver TrackerDriftFMCW',
                current: 0,
                total: 100,
            },
            {
                message: 'Rendering Interference for TrackerDriftFMCW',
                current: 25,
                total: 100,
            },
            {
                message: 'Applying Noise for TrackerDriftFMCW',
                current: 50,
                total: 100,
            },
            {
                message: 'Writing HDF5 for TrackerDriftFMCW',
                current: 75,
                total: 100,
            },
            {
                message: 'Finalized TrackerDriftFMCW',
                current: 100,
                total: 100,
            },
        ]);

        expect(Object.keys(tracked)).toEqual(['receiver:TrackerDriftFMCW']);
        expect(tracked['receiver:TrackerDriftFMCW']).toMatchObject({
            message: 'Finalized TrackerDriftFMCW',
            current: 100,
            total: 100,
            details: [
                {
                    id: 'streaming-finalizing',
                    message: 'Finalizing FMCW Receiver TrackerDriftFMCW',
                },
                {
                    id: 'streaming-rendering-interference',
                    message: 'Rendering Interference for TrackerDriftFMCW',
                },
                {
                    id: 'streaming-applying-noise',
                    message: 'Applying Noise for TrackerDriftFMCW',
                },
                {
                    id: 'streaming-writing-hdf5',
                    message: 'Writing HDF5 for TrackerDriftFMCW',
                },
                {
                    id: 'streaming-finalized',
                    message: 'Finalized TrackerDriftFMCW',
                },
            ],
        });
    });

    test('keeps pulsed export chunks and completion on one stable row', () => {
        const tracked = trackProgress([
            {
                message: 'Exporting PulsedRadar: Chunk 9928',
                current: 9928,
                total: 0,
            },
            {
                message: 'Finished Exporting PulsedRadar',
                current: 100,
                total: 100,
            },
        ]);

        expect(Object.keys(tracked)).toEqual(['receiver:PulsedRadar']);
        expect(tracked['receiver:PulsedRadar']).toEqual({
            message: 'Finished Exporting PulsedRadar',
            current: 100,
            total: 100,
            details: [
                {
                    id: 'export-chunk',
                    message: 'Exporting PulsedRadar: Chunk 9928',
                    current: 9928,
                    total: 0,
                },
                {
                    id: 'export-finished',
                    message: 'Finished Exporting PulsedRadar',
                    current: 100,
                    total: 100,
                },
            ],
        });
    });

    test('keeps only the latest chunk detail for high-volume export updates', () => {
        const tracked = trackProgress([
            {
                message: 'Exporting PulsedRadar: Chunk 10',
                current: 10,
                total: 0,
            },
            {
                message: 'Exporting PulsedRadar: Chunk 9928',
                current: 9928,
                total: 0,
            },
        ]);

        expect(tracked['receiver:PulsedRadar'].details).toEqual([
            {
                id: 'export-chunk',
                message: 'Exporting PulsedRadar: Chunk 9928',
                current: 9928,
                total: 0,
            },
        ]);
    });

    test('keeps main simulation messages out of the detail list', () => {
        const tracked = trackProgress([
            {
                message: 'Simulating... 1.00s / 10.00s',
                current: 10,
                total: 100,
            },
            {
                message: 'Main simulation finished. Waiting for data export...',
                current: 100,
                total: 100,
            },
            {
                message: 'Simulation complete',
                current: 100,
                total: 100,
            },
        ]);

        expect(Object.keys(tracked)).toEqual(['main']);
        expect(tracked.main).toEqual({
            message: 'Simulation complete',
            current: 100,
            total: 100,
        });
    });

    test('marks known stale receiver and export rows complete after success', () => {
        const completed = normalizeCompletedProgressSnapshot({
            CW: {
                message: 'Finalizing CW Receiver CWRadar',
                current: 0,
                total: 100,
            },
            HDF5: {
                message: 'Writing HDF5 for CWRadar',
                current: 75,
                total: 100,
            },
            PulsedRadar: {
                message: 'Exporting PulsedRadar: Chunk 9928',
                current: 9928,
                total: 0,
            },
            main: {
                message: 'Main simulation finished. Waiting for data export...',
                current: 100,
                total: 100,
            },
        });

        expect(completed).toEqual({
            'receiver:CWRadar': {
                message: 'Finalized CWRadar',
                current: 100,
                total: 100,
                details: [
                    {
                        id: 'streaming-finalizing',
                        message: 'Finalizing CW Receiver CWRadar',
                        current: 0,
                        total: 100,
                    },
                    {
                        id: 'streaming-writing-hdf5',
                        message: 'Writing HDF5 for CWRadar',
                        current: 75,
                        total: 100,
                    },
                    {
                        id: 'streaming-finalized',
                        message: 'Finalized CWRadar',
                        current: 100,
                        total: 100,
                    },
                ],
            },
            'receiver:PulsedRadar': {
                message: 'Finished Exporting PulsedRadar',
                current: 100,
                total: 100,
                details: [
                    {
                        id: 'export-chunk',
                        message: 'Exporting PulsedRadar: Chunk 9928',
                        current: 9928,
                        total: 0,
                    },
                    {
                        id: 'export-finished',
                        message: 'Finished Exporting PulsedRadar',
                        current: 100,
                        total: 100,
                    },
                ],
            },
            main: {
                message: 'Simulation complete',
                current: 100,
                total: 100,
            },
        });
    });

    test('clamps determinate progress percentages', () => {
        expect(
            getSimulationProgressPercent({
                message: 'below',
                current: -10,
                total: 100,
            })
        ).toBe(0);
        expect(
            getSimulationProgressPercent({
                message: 'above',
                current: 125,
                total: 100,
            })
        ).toBe(100);
        expect(
            getSimulationProgressPercent({
                message: 'chunk',
                current: 9928,
                total: 0,
            })
        ).toBeNull();
    });
});
