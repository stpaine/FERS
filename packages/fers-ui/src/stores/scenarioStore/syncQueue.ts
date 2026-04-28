// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { invoke } from '@tauri-apps/api/core';

// Single FIFO queue for every backend write, granular or full.
// Serializing all writes through one promise chain prevents the race where
// a granular edit on a freshly-added item would reach the backend before
// the full sync that created it.
let queue: Promise<void> = Promise.resolve();
let invokeBackend: typeof invoke = invoke;

const GRANULAR_FLUSH_INTERVAL_MS = 16;

interface GranularUpdate {
    itemType: string;
    itemId: string;
    json: string;
}

export interface GranularSyncFailure {
    itemType: string;
    itemId: string;
    error: unknown;
}

interface Deferred {
    promise: Promise<void>;
    resolve: () => void;
    reject: (reason?: unknown) => void;
}

type GranularSyncFailureHandler = (
    failure: GranularSyncFailure
) => Promise<void> | void;
type SyncWarningsHandler = (warnings: string[]) => Promise<void> | void;

// Latest pending granular payload per backend object. Older snapshots for the
// same item are obsolete because granular syncs send full item state, not diffs.
const pendingGranularUpdates = new Map<string, GranularUpdate>();
let granularFlushTimer: ReturnType<typeof setTimeout> | null = null;
let scheduledGranularFlush: Deferred | null = null;

// Reference to a full sync that has been enqueued but has not yet started
// executing. While this is set, any further enqueueFullSync() calls are
// coalesced into it — the snapshot is captured when the task runs, so
// later changes are automatically included. It is cleared at the start of
// the task body so subsequent enqueues create a fresh task.
let pendingFullSync: Promise<void> | null = null;
let granularSyncEpoch = 0;
let granularSyncFailureHandler: GranularSyncFailureHandler | null = null;
let syncWarningsHandler: SyncWarningsHandler | null = null;

function createDeferred(): Deferred {
    let resolve!: () => void;
    let reject!: (reason?: unknown) => void;
    const promise = new Promise<void>((res, rej) => {
        resolve = res;
        reject = rej;
    });
    return { promise, resolve, reject };
}

function getGranularUpdateKey(itemType: string, itemId: string): string {
    return `${itemType}:${itemId}`;
}

function discardBufferedGranularSync(replacement?: Promise<void>): void {
    if (granularFlushTimer) {
        clearTimeout(granularFlushTimer);
        granularFlushTimer = null;
    }

    if (pendingGranularUpdates.size > 0) {
        pendingGranularUpdates.clear();
    }

    if (scheduledGranularFlush) {
        const deferred = scheduledGranularFlush;
        if (replacement) {
            void replacement.then(
                () => deferred.resolve(),
                (error) => deferred.reject(error)
            );
        } else {
            deferred.resolve();
        }
        scheduledGranularFlush = null;
    }
}

function invalidatePendingGranularSyncs(): void {
    granularSyncEpoch += 1;
    discardBufferedGranularSync();
}

async function handleGranularSyncFailure(
    update: GranularUpdate,
    error: unknown
): Promise<void> {
    console.error(
        `Granular sync failed (${update.itemType} ${update.itemId}):`,
        error
    );

    invalidatePendingGranularSyncs();

    if (granularSyncFailureHandler) {
        await granularSyncFailureHandler({
            itemType: update.itemType,
            itemId: update.itemId,
            error,
        });
    }
}

async function handleSyncWarnings(warnings: unknown): Promise<void> {
    if (
        !Array.isArray(warnings) ||
        warnings.length === 0 ||
        !syncWarningsHandler
    ) {
        return;
    }
    await syncWarningsHandler(
        warnings.filter(
            (warning): warning is string => typeof warning === 'string'
        )
    );
}

function scheduleGranularFlush(): Promise<void> {
    if (!scheduledGranularFlush) {
        scheduledGranularFlush = createDeferred();
    }

    if (granularFlushTimer) {
        return scheduledGranularFlush.promise;
    }

    granularFlushTimer = setTimeout(() => {
        granularFlushTimer = null;

        const updates = [...pendingGranularUpdates.values()];
        pendingGranularUpdates.clear();

        const flushDeferred = scheduledGranularFlush;
        scheduledGranularFlush = null;

        if (updates.length === 0) {
            flushDeferred?.resolve();
            return;
        }

        const taskEpoch = granularSyncEpoch;
        const task = queue.then(async () => {
            if (taskEpoch !== granularSyncEpoch) {
                return;
            }

            for (const update of updates) {
                try {
                    const warnings = await invokeBackend<string[]>(
                        'update_item_from_json',
                        {
                            itemType: update.itemType,
                            itemId: update.itemId,
                            json: update.json,
                        }
                    );
                    await handleSyncWarnings(warnings);
                } catch (e) {
                    await handleGranularSyncFailure(update, e);
                    throw e;
                }
            }
        });

        queue = task.catch(() => undefined);
        void task.then(
            () => flushDeferred?.resolve(),
            (error) => flushDeferred?.reject(error)
        );
    }, GRANULAR_FLUSH_INTERVAL_MS);

    return scheduledGranularFlush.promise;
}

export function setSyncQueueInvokerForTests(testInvoke?: typeof invoke): void {
    invokeBackend = testInvoke ?? invoke;
}

/** Enqueue a granular item update behind any in-flight work. */
export function enqueueGranularSync(
    itemType: string,
    itemId: string,
    json: string
): Promise<void> {
    if (pendingFullSync) {
        return pendingFullSync;
    }

    pendingGranularUpdates.set(getGranularUpdateKey(itemType, itemId), {
        itemType,
        itemId,
        json,
    });
    return scheduleGranularFlush();
}

export function enqueueGranularSyncDetached(
    itemType: string,
    itemId: string,
    json: string
): void {
    void enqueueGranularSync(itemType, itemId, json).catch(() => undefined);
}

export function registerGranularSyncFailureHandler(
    handler: GranularSyncFailureHandler | null
): void {
    granularSyncFailureHandler = handler;
}

export function registerSyncWarningsHandler(
    handler: SyncWarningsHandler | null
): void {
    syncWarningsHandler = handler;
}

/**
 * Enqueue a full scenario snapshot. Coalesces with any pending full sync.
 * `buildJson` MUST read live state at call time.
 */
export function enqueueFullSync(buildSnapshot: () => string): Promise<void> {
    // A pending full snapshot already contains any granular edits that have not
    // yet been appended to the FIFO queue, so drop those buffered writes.
    if (pendingFullSync) {
        discardBufferedGranularSync(pendingFullSync);
        return pendingFullSync;
    }

    const task = queue.then(async () => {
        // Clear before snapshot so any enqueues from this point create a new task.
        pendingFullSync = null;
        const json = buildSnapshot();
        try {
            const warnings = await invokeBackend<string[]>(
                'update_scenario_from_json',
                { json }
            );
            await handleSyncWarnings(warnings);
        } catch (e) {
            console.error('Full sync failed:', e);
        }
    });
    pendingFullSync = task;
    discardBufferedGranularSync(task);
    queue = task.catch(() => undefined);
    return task;
}

/** Resolves once every currently-queued sync task has settled. */
export function waitForSyncIdle(): Promise<void> {
    return Promise.all([
        queue,
        scheduledGranularFlush?.promise ?? Promise.resolve(),
    ]).then(() => undefined);
}

export function resetSyncQueueForTests(): void {
    if (granularFlushTimer) {
        clearTimeout(granularFlushTimer);
        granularFlushTimer = null;
    }

    pendingGranularUpdates.clear();
    scheduledGranularFlush = null;
    pendingFullSync = null;
    queue = Promise.resolve();
    invokeBackend = invoke;
    granularSyncEpoch = 0;
    granularSyncFailureHandler = null;
    syncWarningsHandler = null;
}
