// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { invoke } from '@tauri-apps/api/core';

// Single FIFO queue for every backend write, granular or full.
// Serializing all writes through one promise chain prevents the race where
// a granular edit on a freshly-added item would reach the backend before
// the full sync that created it.
let queue: Promise<void> = Promise.resolve();

// Reference to a full sync that has been enqueued but has not yet started
// executing. While this is set, any further enqueueFullSync() calls are
// coalesced into it — the snapshot is captured when the task runs, so
// later changes are automatically included. It is cleared at the start of
// the task body so subsequent enqueues create a fresh task.
let pendingFullSync: Promise<void> | null = null;

/** Enqueue a granular item update behind any in-flight work. */
export function enqueueGranularSync(
    itemType: string,
    itemId: string,
    json: string
): Promise<void> {
    const task = queue.then(async () => {
        try {
            await invoke('update_item_from_json', { itemType, itemId, json });
        } catch (e) {
            console.error(`Granular sync failed (${itemType} ${itemId}):`, e);
            throw e;
        }
    });
    queue = task.catch(() => undefined);
    return task;
}

/**
 * Enqueue a full scenario snapshot. Coalesces with any pending full sync.
 * `buildJson` MUST read live state at call time.
 */
export function enqueueFullSync(buildSnapshot: () => string): Promise<void> {
    if (pendingFullSync) return pendingFullSync;

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
    queue = task.catch(() => undefined);
    return task;
}

/** Resolves once every currently-queued sync task has settled. */
export function waitForSyncIdle(): Promise<void> {
    return queue;
}
