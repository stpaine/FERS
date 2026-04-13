// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { create } from 'zustand';

export type SimulationProgressState = {
    message: string;
    current: number;
    total: number;
};

export type SimulationRunStatus = 'idle' | 'running' | 'completed' | 'failed';

type SimulationProgressStore = {
    isSimulating: boolean;
    isGeneratingKml: boolean;
    simulationProgress: Record<string, SimulationProgressState>;
    simulationRunStatus: SimulationRunStatus;
    simulationRunError: string | null;

    setIsSimulating: (isSimulating: boolean) => void;
    setIsGeneratingKml: (isGeneratingKml: boolean) => void;
    startSimulationRun: () => void;
    setSimulationProgressSnapshot: (
        progress: Record<string, SimulationProgressState>
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

        setIsSimulating: (isSimulating) => set({ isSimulating }),
        setIsGeneratingKml: (isGeneratingKml) => set({ isGeneratingKml }),
        startSimulationRun: () =>
            set({
                isSimulating: true,
                simulationProgress: {},
                simulationRunStatus: 'running',
                simulationRunError: null,
            }),
        setSimulationProgressSnapshot: (progress) =>
            set({ simulationProgress: progress }),
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
            }),
    })
);
