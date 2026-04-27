// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import type {
    SimulationProgressDetail,
    SimulationProgressState,
} from '@/stores/simulationProgressStore';

type ProgressMessageMatch =
    | { kind: 'main'; key: 'main'; detailId: string }
    | { kind: 'streaming'; key: string; receiverName: string; detailId: string }
    | { kind: 'export'; key: string; receiverName: string; detailId: string }
    | { kind: 'unknown'; key: string; detailId: string };

export type NormalizedSimulationProgress = {
    key: string;
    progress: SimulationProgressState;
};

const receiverKey = (receiverName: string) => `receiver:${receiverName}`;

const detailOrder: Record<string, number> = {
    'main-progress': 0,
    'main-export-wait': 1,
    'main-complete': 2,
    'streaming-finalizing': 0,
    'streaming-rendering-interference': 1,
    'streaming-applying-noise': 2,
    'streaming-writing-hdf5': 3,
    'streaming-finalized': 4,
    'export-chunk': 0,
    'export-finished': 1,
};

const sortDetails = (details: SimulationProgressDetail[]) =>
    [...details].sort(
        (a, b) =>
            (detailOrder[a.id] ?? Number.MAX_SAFE_INTEGER) -
            (detailOrder[b.id] ?? Number.MAX_SAFE_INTEGER)
    );

const buildReceiverMatch = (
    kind: 'streaming' | 'export',
    receiverName: string,
    detailId: string
): ProgressMessageMatch => ({
    kind,
    key: receiverKey(receiverName),
    receiverName,
    detailId,
});

const stripDetails = (
    progress: SimulationProgressState | SimulationProgressDetail
): SimulationProgressState => ({
    message: progress.message,
    current: progress.current,
    total: progress.total,
});

const matchProgressMessage = (message: string): ProgressMessageMatch => {
    if (
        message.startsWith('Simulating') ||
        message.startsWith('Initializing')
    ) {
        return { kind: 'main', key: 'main', detailId: 'main-progress' };
    }

    if (message.startsWith('Main simulation finished')) {
        return { kind: 'main', key: 'main', detailId: 'main-export-wait' };
    }

    if (message === 'Simulation complete') {
        return { kind: 'main', key: 'main', detailId: 'main-complete' };
    }

    if (message.startsWith('Finalizing CW Receiver ')) {
        return buildReceiverMatch(
            'streaming',
            message.slice('Finalizing CW Receiver '.length),
            'streaming-finalizing'
        );
    }

    if (message.startsWith('Finalizing FMCW Receiver ')) {
        return buildReceiverMatch(
            'streaming',
            message.slice('Finalizing FMCW Receiver '.length),
            'streaming-finalizing'
        );
    }

    if (message.startsWith('Rendering Interference for ')) {
        return buildReceiverMatch(
            'streaming',
            message.slice('Rendering Interference for '.length),
            'streaming-rendering-interference'
        );
    }

    if (message.startsWith('Applying Noise for ')) {
        return buildReceiverMatch(
            'streaming',
            message.slice('Applying Noise for '.length),
            'streaming-applying-noise'
        );
    }

    if (message.startsWith('Writing HDF5 for ')) {
        return buildReceiverMatch(
            'streaming',
            message.slice('Writing HDF5 for '.length),
            'streaming-writing-hdf5'
        );
    }

    if (message.startsWith('Finalized ')) {
        return buildReceiverMatch(
            'streaming',
            message.slice('Finalized '.length),
            'streaming-finalized'
        );
    }

    const exportChunkMatch = message.match(/^Exporting (.+): Chunk \d+$/);
    if (exportChunkMatch) {
        return buildReceiverMatch(
            'export',
            exportChunkMatch[1],
            'export-chunk'
        );
    }

    if (message.startsWith('Finished Exporting ')) {
        return buildReceiverMatch(
            'export',
            message.slice('Finished Exporting '.length),
            'export-finished'
        );
    }

    return { kind: 'unknown', key: message, detailId: message };
};

export const normalizeSimulationProgressEvent = (
    progress: SimulationProgressState
): NormalizedSimulationProgress => {
    const match = matchProgressMessage(progress.message);

    return {
        key: match.key,
        progress,
    };
};

const toDetail = (
    match: ProgressMessageMatch,
    progress: SimulationProgressState
): SimulationProgressDetail => ({
    id: match.detailId,
    ...stripDetails(progress),
});

const upsertDetail = (
    details: SimulationProgressDetail[],
    detail: SimulationProgressDetail
) => {
    const detailIndex = details.findIndex((item) => item.id === detail.id);

    if (detailIndex === -1) {
        return sortDetails([...details, detail]);
    }

    return sortDetails(
        details.map((item, index) => (index === detailIndex ? detail : item))
    );
};

const mergeProgress = (
    existing: SimulationProgressState | undefined,
    progress: SimulationProgressState
): SimulationProgressState => {
    const match = matchProgressMessage(progress.message);
    const existingDetails = existing?.details ?? [];

    if (match.kind === 'main') {
        return stripDetails(progress);
    }

    return {
        ...stripDetails(progress),
        details: upsertDetail(existingDetails, toDetail(match, progress)),
    };
};

export const addSimulationProgressEvent = (
    progress: Record<string, SimulationProgressState>,
    event: SimulationProgressState
): Record<string, SimulationProgressState> => {
    const { key, progress: normalizedProgress } =
        normalizeSimulationProgressEvent(event);

    return {
        ...progress,
        [key]: mergeProgress(progress[key], normalizedProgress),
    };
};

const completeProgressForMatch = (
    match: ProgressMessageMatch,
    progress: SimulationProgressState
): SimulationProgressState => {
    if (match.kind === 'main') {
        return {
            message: 'Simulation complete',
            current: 100,
            total: 100,
        };
    }

    if (match.kind === 'streaming') {
        return {
            message: `Finalized ${match.receiverName}`,
            current: 100,
            total: 100,
        };
    }

    if (match.kind === 'export') {
        return {
            message: `Finished Exporting ${match.receiverName}`,
            current: 100,
            total: 100,
        };
    }

    return progress;
};

const getProgressSortValue = (progress: SimulationProgressState) =>
    getSimulationProgressPercent(progress) ?? 0;

const chooseLatestProgress = (
    existing: SimulationProgressState | undefined,
    next: SimulationProgressState
) => {
    if (!existing) {
        return next;
    }

    return getProgressSortValue(next) >= getProgressSortValue(existing)
        ? next
        : existing;
};

export const normalizeCompletedProgressSnapshot = (
    progress: Record<string, SimulationProgressState>
): Record<string, SimulationProgressState> => {
    const completedProgress: Record<string, SimulationProgressState> = {};

    for (const item of Object.values(progress)) {
        for (const detail of item.details ?? []) {
            const detailProgress = stripDetails(detail);
            const { key } = normalizeSimulationProgressEvent(detailProgress);
            completedProgress[key] = mergeProgress(
                completedProgress[key],
                detailProgress
            );
        }

        const originalItem = stripDetails(item);
        const match = matchProgressMessage(item.message);
        completedProgress[match.key] = mergeProgress(
            completedProgress[match.key],
            originalItem
        );

        const completedItem = completeProgressForMatch(match, originalItem);
        const mergedItem = mergeProgress(
            completedProgress[match.key],
            completedItem
        );
        completedProgress[match.key] = chooseLatestProgress(
            completedProgress[match.key],
            mergedItem
        );
    }

    completedProgress.main = mergeProgress(completedProgress.main, {
        message: 'Simulation complete',
        current: 100,
        total: 100,
    });

    return completedProgress;
};

export const getSimulationProgressPercent = (
    progress: SimulationProgressState
): number | null => {
    if (
        progress.total <= 0 ||
        !Number.isFinite(progress.current) ||
        !Number.isFinite(progress.total)
    ) {
        return null;
    }

    const percent = (progress.current / progress.total) * 100;

    return Math.min(100, Math.max(0, percent));
};
