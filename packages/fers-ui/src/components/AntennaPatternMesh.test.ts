// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { describe, expect, it } from 'bun:test';
import {
    BACKEND_BUSY_MESSAGE,
    getAntennaPreviewErrorAction,
    isBackendBusyError,
    shouldDeferAntennaPreviewFetch,
} from './AntennaPatternMesh';

describe('AntennaPatternMesh preview recovery', () => {
    it('treats backend lock contention as a transient retry condition', () => {
        const error = new Error(BACKEND_BUSY_MESSAGE);

        expect(isBackendBusyError(error)).toBe(true);
        expect(getAntennaPreviewErrorAction(error)).toEqual({
            clearError: true,
            clearPattern: false,
            scheduleRetry: true,
        });
    });

    it('preserves hard failures as visible preview errors', () => {
        const error = new Error('antenna pattern file missing');

        expect(isBackendBusyError(error)).toBe(false);
        expect(getAntennaPreviewErrorAction(error)).toEqual({
            clearError: false,
            clearPattern: true,
            scheduleRetry: false,
        });
    });

    it('defers fetches until backend work has finished', () => {
        expect(
            shouldDeferAntennaPreviewFetch({
                hasRequest: true,
                isBackendSyncing: true,
                isSimulating: false,
                isGeneratingKml: false,
            })
        ).toBe(true);

        expect(
            shouldDeferAntennaPreviewFetch({
                hasRequest: true,
                isBackendSyncing: false,
                isSimulating: true,
                isGeneratingKml: false,
            })
        ).toBe(true);

        expect(
            shouldDeferAntennaPreviewFetch({
                hasRequest: true,
                isBackendSyncing: false,
                isSimulating: false,
                isGeneratingKml: true,
            })
        ).toBe(true);

        expect(
            shouldDeferAntennaPreviewFetch({
                hasRequest: true,
                isBackendSyncing: false,
                isSimulating: false,
                isGeneratingKml: false,
            })
        ).toBe(false);

        expect(
            shouldDeferAntennaPreviewFetch({
                hasRequest: false,
                isBackendSyncing: true,
                isSimulating: true,
                isGeneratingKml: true,
            })
        ).toBe(false);
    });
});
