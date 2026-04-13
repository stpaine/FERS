// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { beforeEach, describe, expect, test } from 'bun:test';
import { useSimulationProgressStore } from '../simulationProgressStore';
import { useScenarioStore } from './index';
import {
    resetSyncQueueForTests,
    setSyncQueueInvokerForTests,
    waitForSyncIdle,
} from './syncQueue';

type InvokeFn = typeof import('@tauri-apps/api/core').invoke;

describe('simulation progress lifecycle', () => {
    beforeEach(() => {
        resetSyncQueueForTests();
        setSyncQueueInvokerForTests((async () => []) as InvokeFn);
        useSimulationProgressStore.setState({
            isSimulating: false,
            simulationProgress: {},
            simulationRunStatus: 'idle',
            simulationRunError: null,
            simulationOutputMetadata: null,
        });
    });

    test('clears previous progress when a simulation starts', () => {
        useSimulationProgressStore.setState({
            simulationProgress: {
                main: {
                    message: 'Finished previous simulation',
                    current: 100,
                    total: 100,
                },
            },
            simulationRunStatus: 'completed',
            simulationRunError: 'old error',
            simulationOutputMetadata: {
                schema_version: 1,
                simulation_name: 'previous',
                output_directory: '.',
                start_time: 0,
                end_time: 1,
                sampling_rate: 1000,
                oversample_ratio: 1,
                files: [],
            },
        });

        useSimulationProgressStore.getState().startSimulationRun();

        const state = useSimulationProgressStore.getState();
        expect(state.isSimulating).toBe(true);
        expect(state.simulationProgress).toEqual({});
        expect(state.simulationRunStatus).toBe('running');
        expect(state.simulationRunError).toBeNull();
        expect(state.simulationOutputMetadata).toBeNull();
    });

    test('preserves progress when a simulation completes', () => {
        const progress = {
            main: {
                message: 'Simulating: 100%',
                current: 100,
                total: 100,
            },
            Receiver1: {
                message: 'Finished Exporting Receiver1',
                current: 100,
                total: 100,
            },
        };

        useSimulationProgressStore.getState().startSimulationRun();
        useSimulationProgressStore
            .getState()
            .setSimulationProgressSnapshot(progress);
        useSimulationProgressStore.getState().completeSimulationRun();

        const state = useSimulationProgressStore.getState();
        expect(state.isSimulating).toBe(false);
        expect(state.simulationRunStatus).toBe('completed');
        expect(state.simulationRunError).toBeNull();
        expect(state.simulationProgress).toEqual(progress);
    });

    test('preserves progress and records the error when a simulation fails', () => {
        const progress = {
            main: {
                message: 'Simulating: 40%',
                current: 40,
                total: 100,
            },
        };

        useSimulationProgressStore.getState().startSimulationRun();
        useSimulationProgressStore
            .getState()
            .setSimulationProgressSnapshot(progress);
        useSimulationProgressStore
            .getState()
            .failSimulationRun('Simulation failed');

        const state = useSimulationProgressStore.getState();
        expect(state.isSimulating).toBe(false);
        expect(state.simulationRunStatus).toBe('failed');
        expect(state.simulationRunError).toBe('Simulation failed');
        expect(state.simulationProgress).toEqual(progress);
    });

    test('resetScenario clears stored simulation progress', async () => {
        useSimulationProgressStore.setState({
            isSimulating: true,
            simulationProgress: {
                main: {
                    message: 'Simulating: 70%',
                    current: 70,
                    total: 100,
                },
            },
            simulationRunStatus: 'running',
            simulationRunError: 'old error',
            simulationOutputMetadata: {
                schema_version: 1,
                simulation_name: 'old',
                output_directory: '.',
                start_time: 0,
                end_time: 1,
                sampling_rate: 1000,
                oversample_ratio: 1,
                files: [],
            },
        });

        useScenarioStore.getState().resetScenario();
        await waitForSyncIdle();

        const state = useSimulationProgressStore.getState();
        expect(state.isSimulating).toBe(false);
        expect(state.simulationProgress).toEqual({});
        expect(state.simulationRunStatus).toBe('idle');
        expect(state.simulationRunError).toBeNull();
        expect(state.simulationOutputMetadata).toBeNull();
    });
});
