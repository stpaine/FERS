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
    | 'OFF'
    | 'UNKNOWN';

export type FersLogEntry = {
    sequence: number;
    level: FersLogLevel;
    line: string;
};

export type ConfigurableFersLogLevel = Exclude<FersLogLevel, 'UNKNOWN'>;

type FersLogStore = {
    entries: FersLogEntry[];
    droppedCount: number;
    maxLines: number;
    logLevel: ConfigurableFersLogLevel;
    isOpen: boolean;
    appendLog: (entry: FersLogEntry) => void;
    clearLogs: () => void;
    setMaxLines: (maxLines: number) => void;
    setLogLevel: (logLevel: ConfigurableFersLogLevel) => void;
    setOpen: (isOpen: boolean) => void;
    toggleOpen: () => void;
};

export const DEFAULT_LOG_MAX_LINES = 2000;
export const MIN_LOG_MAX_LINES = 100;
export const MAX_LOG_MAX_LINES = 20_000;
export const DEFAULT_LOG_LEVEL: ConfigurableFersLogLevel = 'INFO';
export const LOG_LEVEL_OPTIONS: ConfigurableFersLogLevel[] = [
    'TRACE',
    'DEBUG',
    'INFO',
    'WARNING',
    'ERROR',
    'FATAL',
    'OFF',
];

export const isConfigurableFersLogLevel = (
    level: FersLogLevel
): level is ConfigurableFersLogLevel =>
    LOG_LEVEL_OPTIONS.some((option) => option === level);

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
    logLevel: DEFAULT_LOG_LEVEL,
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
    setLogLevel: (logLevel) => set({ logLevel }),
    setOpen: (isOpen) => set({ isOpen }),
    toggleOpen: () => set((state) => ({ isOpen: !state.isOpen })),
}));
