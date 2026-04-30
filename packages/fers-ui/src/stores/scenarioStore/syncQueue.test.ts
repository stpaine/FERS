// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { beforeEach, describe, expect, test } from 'bun:test';
import {
    enqueueFullSync,
    enqueueGranularSync,
    type GranularSyncFailure,
    registerGranularSyncFailureHandler,
    resetSyncQueueForTests,
    setSyncQueueInvokerForTests,
    waitForSyncIdle,
} from './syncQueue';

type InvokeFn = typeof import('@tauri-apps/api/core').invoke;

function delay(ms: number): Promise<void> {
    return new Promise((resolve) => setTimeout(resolve, ms));
}

function createDeferred() {
    let resolve!: () => void;
    const promise = new Promise<void>((res) => {
        resolve = res;
    });
    return { promise, resolve };
}

describe('syncQueue granular recovery', () => {
    beforeEach(() => {
        resetSyncQueueForTests();
    });

    test('does not invoke recovery when granular sync succeeds', async () => {
        const invocations: string[] = [];
        const failures: GranularSyncFailure[] = [];

        setSyncQueueInvokerForTests((async (
            command: string,
            args?: Record<string, unknown>
        ) => {
            invocations.push(`${command}:${String(args?.itemId ?? '')}`);
        }) as InvokeFn);
        registerGranularSyncFailureHandler((failure) => {
            failures.push(failure);
        });

        await enqueueGranularSync('Waveform', '10', '{"id":"10"}');
        await waitForSyncIdle();

        expect(invocations).toEqual(['update_item_from_json:10']);
        expect(failures).toHaveLength(0);
    });

    test('coalesces granular edits into a pending full sync', async () => {
        const invocations: string[] = [];
        let snapshot = '{"simulation":{"name":"before"}}';

        setSyncQueueInvokerForTests((async (
            command: string,
            args?: Record<string, unknown>
        ) => {
            invocations.push(
                `${command}:${String(args?.itemId ?? '')}:${String(args?.json ?? '')}`
            );
        }) as InvokeFn);

        const fullSync = enqueueFullSync(() => snapshot);
        snapshot = '{"simulation":{"name":"after"}}';
        const granularSync = enqueueGranularSync(
            'Platform',
            '281474976710657',
            '{"id":"281474976710657","name":"after"}'
        );

        await fullSync;
        await granularSync;
        await waitForSyncIdle();

        expect(invocations).toEqual([
            'update_scenario_from_json::{"simulation":{"name":"after"}}',
        ]);
    });

    test('rejects failed full syncs', async () => {
        setSyncQueueInvokerForTests((async (command: string) => {
            if (command === 'update_scenario_from_json') {
                throw new Error('backend rejected scenario');
            }
        }) as InvokeFn);

        await expect(
            enqueueFullSync(() => '{"simulation":{"name":"bad"}}')
        ).rejects.toThrow('backend rejected scenario');
        await waitForSyncIdle();
    });

    test('recovers once and discards stale queued granular flushes after a failure', async () => {
        const invocations: string[] = [];
        const failures: GranularSyncFailure[] = [];
        const firstFailureGate = createDeferred();
        const recoveryGate = createDeferred();
        let signalRecoveryStarted!: () => void;
        const recoveryStarted = new Promise<void>((resolve) => {
            signalRecoveryStarted = resolve;
        });

        setSyncQueueInvokerForTests((async (
            command: string,
            args?: Record<string, unknown>
        ) => {
            const itemId = String(args?.itemId ?? '');
            invocations.push(`${command}:${itemId}`);

            if (command === 'update_item_from_json' && itemId === '1') {
                await firstFailureGate.promise;
                throw new Error('backend rejected update');
            }
        }) as InvokeFn);
        registerGranularSyncFailureHandler(async (failure) => {
            failures.push(failure);
            signalRecoveryStarted();
            await recoveryGate.promise;
        });

        const failingFlush = enqueueGranularSync('Waveform', '1', '{"id":"1"}');
        await delay(25);

        const staleFlush = enqueueGranularSync('Waveform', '2', '{"id":"2"}');
        await delay(25);

        firstFailureGate.resolve();
        await recoveryStarted;

        const freshFlush = enqueueGranularSync('Waveform', '3', '{"id":"3"}');
        await delay(25);

        recoveryGate.resolve();

        await expect(failingFlush).rejects.toThrow('backend rejected update');
        await staleFlush;
        await freshFlush;
        await waitForSyncIdle();

        expect(failures).toHaveLength(1);
        expect(failures[0]?.itemId).toBe('1');
        expect(invocations).toEqual([
            'update_item_from_json:1',
            'update_item_from_json:3',
        ]);
    });
});
