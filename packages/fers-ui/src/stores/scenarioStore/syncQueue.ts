// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { invoke } from '@tauri-apps/api/core';

// Single FIFO queue for every backend write, granular or full.
// Serializing all writes through one promise chain prevents the race where
// a granular edit on a freshly-added item would reach the backend before
// the full sync that created it.
let queue: Promise<void> = Promise.resolve();

const GRANULAR_FLUSH_INTERVAL_MS = 16;

interface GranularUpdate {
    itemType: string;
    itemId: string;
    json: string;
}

interface Deferred {
    promise: Promise<void>;
    resolve: () => void;
    reject: (reason?: unknown) => void;
}

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

        const task = queue.then(async () => {
            let firstError: unknown = null;

            for (const update of updates) {
                try {
                    await invoke('update_item_from_json', {
                        itemType: update.itemType,
                        itemId: update.itemId,
                        json: update.json,
                    });
                } catch (e) {
                    console.error(
                        `Granular sync failed (${update.itemType} ${update.itemId}):`,
                        e
                    );
                    firstError ??= e;
                }
            }

            if (firstError) {
                throw firstError;
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

/** Enqueue a granular item update behind any in-flight work. */
export function enqueueGranularSync(
    itemType: string,
    itemId: string,
    json: string
): Promise<void> {
    pendingGranularUpdates.set(getGranularUpdateKey(itemType, itemId), {
        itemType,
        itemId,
        json,
    });
    return scheduleGranularFlush();
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
            await invoke('update_scenario_from_json', { json });
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
