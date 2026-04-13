// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { create } from 'zustand';

export type FersLogLevel =
    | 'TRACE'
    | 'DEBUG'
    | 'INFO'
    | 'WARNING'
    | 'ERROR'
    | 'FATAL'
    | 'UNKNOWN';

export type FersLogEntry = {
    sequence: number;
    level: FersLogLevel;
    line: string;
};

type FersLogStore = {
    entries: FersLogEntry[];
    droppedCount: number;
    maxLines: number;
    isOpen: boolean;
    appendLog: (entry: FersLogEntry) => void;
    clearLogs: () => void;
    setMaxLines: (maxLines: number) => void;
    setOpen: (isOpen: boolean) => void;
    toggleOpen: () => void;
};

export const DEFAULT_LOG_MAX_LINES = 2000;
export const MIN_LOG_MAX_LINES = 100;
export const MAX_LOG_MAX_LINES = 20_000;

export const clampLogMaxLines = (maxLines: number) => {
    if (!Number.isFinite(maxLines)) {
        return DEFAULT_LOG_MAX_LINES;
    }

    return Math.max(
        MIN_LOG_MAX_LINES,
        Math.min(MAX_LOG_MAX_LINES, Math.floor(maxLines))
    );
};

const trimEntries = (
    entries: FersLogEntry[],
    maxLines: number
): { entries: FersLogEntry[]; dropped: number } => {
    if (entries.length <= maxLines) {
        return { entries, dropped: 0 };
    }

    return {
        entries: entries.slice(entries.length - maxLines),
        dropped: entries.length - maxLines,
    };
};

export const useFersLogStore = create<FersLogStore>()((set) => ({
    entries: [],
    droppedCount: 0,
    maxLines: DEFAULT_LOG_MAX_LINES,
    isOpen: false,
    appendLog: (entry) =>
        set((state) => {
            const trimmed = trimEntries(
                [...state.entries, entry],
                state.maxLines
            );

            return {
                entries: trimmed.entries,
                droppedCount: state.droppedCount + trimmed.dropped,
            };
        }),
    clearLogs: () => set({ entries: [], droppedCount: 0 }),
    setMaxLines: (maxLines) =>
        set((state) => {
            const nextMaxLines = clampLogMaxLines(maxLines);
            const trimmed = trimEntries(state.entries, nextMaxLines);

            return {
                maxLines: nextMaxLines,
                entries: trimmed.entries,
                droppedCount: state.droppedCount + trimmed.dropped,
            };
        }),
    setOpen: (isOpen) => set({ isOpen }),
    toggleOpen: () => set((state) => ({ isOpen: !state.isOpen })),
}));
