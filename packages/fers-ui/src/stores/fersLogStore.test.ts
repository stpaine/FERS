// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import {
    afterAll,
    beforeAll,
    beforeEach,
    describe,
    expect,
    test,
} from 'bun:test';
import {
    clampLogDrawerWidth,
    clampLogMaxLines,
    DEFAULT_LOG_DRAWER_WIDTH,
    DEFAULT_LOG_LEVEL,
    DEFAULT_LOG_MAX_LINES,
    getEffectiveLogDrawerWidth,
    MAX_LOG_DRAWER_WIDTH,
    MAX_LOG_MAX_LINES,
    MIN_LOG_DRAWER_WIDTH,
    MIN_LOG_MAX_LINES,
    readStoredLogDrawerWidth,
    useFersLogStore,
} from './fersLogStore';

type StorageMock = {
    clear: () => void;
    getItem: (key: string) => string | null;
    removeItem: (key: string) => void;
    setItem: (key: string, value: string) => void;
};

function createStorageMock(): StorageMock {
    const values = new Map<string, string>();

    return {
        clear: () => values.clear(),
        getItem: (key) => values.get(key) ?? null,
        removeItem: (key) => values.delete(key),
        setItem: (key, value) => {
            values.set(key, value);
        },
    };
}

const localStorageMock = createStorageMock();
const originalLocalStorage = globalThis.localStorage;

describe('fers log store', () => {
    beforeAll(() => {
        Object.defineProperty(globalThis, 'localStorage', {
            value: localStorageMock,
            configurable: true,
        });
    });

    afterAll(() => {
        Object.defineProperty(globalThis, 'localStorage', {
            value: originalLocalStorage,
            configurable: true,
        });
    });

    beforeEach(() => {
        globalThis.localStorage.clear();
        useFersLogStore.setState({
            entries: [],
            droppedCount: 0,
            maxLines: DEFAULT_LOG_MAX_LINES,
            logLevel: DEFAULT_LOG_LEVEL,
            isOpen: false,
            drawerWidth: DEFAULT_LOG_DRAWER_WIDTH,
        });
    });

    test('appends log entries in order', () => {
        useFersLogStore
            .getState()
            .appendLog({ sequence: 1, level: 'INFO', line: 'first' });
        useFersLogStore
            .getState()
            .appendLog({ sequence: 2, level: 'ERROR', line: 'second' });

        expect(useFersLogStore.getState().entries).toEqual([
            { sequence: 1, level: 'INFO', line: 'first' },
            { sequence: 2, level: 'ERROR', line: 'second' },
        ]);
    });

    test('trims to bounded buffer and records dropped lines', () => {
        useFersLogStore.setState({ maxLines: 2 });

        for (let sequence = 1; sequence <= 4; sequence += 1) {
            useFersLogStore.getState().appendLog({
                sequence,
                level: 'INFO',
                line: `line ${sequence}`,
            });
        }

        const state = useFersLogStore.getState();
        expect(state.entries.map((entry) => entry.sequence)).toEqual([3, 4]);
        expect(state.droppedCount).toBe(2);
    });

    test('clear resets entries and dropped count', () => {
        useFersLogStore.setState({
            entries: [{ sequence: 1, level: 'INFO', line: 'line' }],
            droppedCount: 4,
        });

        useFersLogStore.getState().clearLogs();

        expect(useFersLogStore.getState().entries).toEqual([]);
        expect(useFersLogStore.getState().droppedCount).toBe(0);
    });

    test('max lines setting clamps and trims immediately', () => {
        useFersLogStore.setState({
            entries: [
                { sequence: 1, level: 'INFO', line: 'one' },
                { sequence: 2, level: 'INFO', line: 'two' },
            ],
            maxLines: 2,
        });

        useFersLogStore.getState().setMaxLines(1);

        const state = useFersLogStore.getState();
        expect(state.maxLines).toBe(MIN_LOG_MAX_LINES);
        expect(state.entries.map((entry) => entry.sequence)).toEqual([1, 2]);
    });

    test('clampLogMaxLines handles invalid and out-of-range values', () => {
        expect(clampLogMaxLines(Number.NaN)).toBe(DEFAULT_LOG_MAX_LINES);
        expect(clampLogMaxLines(1)).toBe(MIN_LOG_MAX_LINES);
        expect(clampLogMaxLines(MAX_LOG_MAX_LINES + 1)).toBe(MAX_LOG_MAX_LINES);
        expect(clampLogMaxLines(1234.9)).toBe(1234);
    });

    test('clampLogDrawerWidth handles invalid and out-of-range values', () => {
        expect(clampLogDrawerWidth(Number.NaN)).toBe(DEFAULT_LOG_DRAWER_WIDTH);
        expect(clampLogDrawerWidth(1)).toBe(MIN_LOG_DRAWER_WIDTH);
        expect(clampLogDrawerWidth(MAX_LOG_DRAWER_WIDTH + 1)).toBe(
            MAX_LOG_DRAWER_WIDTH
        );
        expect(clampLogDrawerWidth(733.6)).toBe(734);
    });

    test('log level defaults and updates', () => {
        expect(useFersLogStore.getState().logLevel).toBe(DEFAULT_LOG_LEVEL);

        useFersLogStore.getState().setLogLevel('TRACE');
        expect(useFersLogStore.getState().logLevel).toBe('TRACE');

        useFersLogStore.getState().setLogLevel('OFF');
        expect(useFersLogStore.getState().logLevel).toBe('OFF');
    });

    test('setDrawerWidth updates state and persists width', () => {
        useFersLogStore.getState().setDrawerWidth(734);

        expect(useFersLogStore.getState().drawerWidth).toBe(734);
        expect(
            JSON.parse(
                globalThis.localStorage.getItem('fers.logViewer.v1') ?? '{}'
            )
        ).toEqual({ drawerWidth: 734 });
    });

    test('readStoredLogDrawerWidth hydrates persisted width', () => {
        globalThis.localStorage.setItem(
            'fers.logViewer.v1',
            JSON.stringify({ drawerWidth: 1440 })
        );

        expect(readStoredLogDrawerWidth()).toBe(MAX_LOG_DRAWER_WIDTH);
    });

    test('effective width clamps to viewport without mutating preferred width', () => {
        useFersLogStore.getState().setDrawerWidth(800);

        const preferredWidth = useFersLogStore.getState().drawerWidth;
        const effectiveWidth = getEffectiveLogDrawerWidth(preferredWidth, 500);

        expect(preferredWidth).toBe(800);
        expect(effectiveWidth).toBe(500);
    });
});
