// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { beforeEach, describe, expect, test } from 'bun:test';
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
        useScenarioStore.setState({
            isSimulating: false,
            simulationProgress: {},
            simulationRunStatus: 'idle',
            simulationRunError: null,
        });
    });

    test('clears previous progress when a simulation starts', () => {
        useScenarioStore.setState({
            simulationProgress: {
                main: {
                    message: 'Finished previous simulation',
                    current: 100,
                    total: 100,
                },
            },
            simulationRunStatus: 'completed',
            simulationRunError: 'old error',
        });

        useScenarioStore.getState().startSimulationRun();

        const state = useScenarioStore.getState();
        expect(state.isSimulating).toBe(true);
        expect(state.simulationProgress).toEqual({});
        expect(state.simulationRunStatus).toBe('running');
        expect(state.simulationRunError).toBeNull();
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

        useScenarioStore.getState().startSimulationRun();
        useScenarioStore.getState().setSimulationProgressSnapshot(progress);
        useScenarioStore.getState().completeSimulationRun();

        const state = useScenarioStore.getState();
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

        useScenarioStore.getState().startSimulationRun();
        useScenarioStore.getState().setSimulationProgressSnapshot(progress);
        useScenarioStore.getState().failSimulationRun('Simulation failed');

        const state = useScenarioStore.getState();
        expect(state.isSimulating).toBe(false);
        expect(state.simulationRunStatus).toBe('failed');
        expect(state.simulationRunError).toBe('Simulation failed');
        expect(state.simulationProgress).toEqual(progress);
    });

    test('resetScenario clears stored simulation progress', async () => {
        useScenarioStore.setState({
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
        });

        useScenarioStore.getState().resetScenario();
        await waitForSyncIdle();

        const state = useScenarioStore.getState();
        expect(state.isSimulating).toBe(false);
        expect(state.simulationProgress).toEqual({});
        expect(state.simulationRunStatus).toBe('idle');
        expect(state.simulationRunError).toBeNull();
    });
});
