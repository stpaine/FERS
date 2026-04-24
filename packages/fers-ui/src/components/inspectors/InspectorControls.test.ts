// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import {
    formatNumberFieldValue,
    resolveNumberFieldBlur,
    resolveTextFieldBlur,
} from './InspectorControls';

describe('InspectorControls blur resolution', () => {
    test('reverts required numeric fields when blurred empty', () => {
        expect(resolveNumberFieldBlur('', 42, 'revert')).toEqual({
            action: 'revert',
            error: 'Value required.',
            nextDraft: '42',
        });
    });

    test('commits null for unsettable numeric fields when blurred empty', () => {
        expect(resolveNumberFieldBlur('', 42, 'null')).toEqual({
            action: 'commit',
            nextDraft: '',
            value: null,
        });
    });

    test('reverts timing noise-entry drafts instead of producing empty commits', () => {
        expect(resolveNumberFieldBlur('', 0.5, 'revert')).toEqual({
            action: 'revert',
            error: 'Value required.',
            nextDraft: '0.5',
        });
    });

    test('rejects invalid numeric drafts on blur', () => {
        expect(resolveNumberFieldBlur('-', 7, 'revert')).toEqual({
            action: 'revert',
            error: 'Invalid number.',
            nextDraft: '7',
        });
    });

    test('parses valid numeric drafts on blur', () => {
        expect(resolveNumberFieldBlur('1e3', 7, 'revert')).toEqual({
            action: 'commit',
            nextDraft: '1000',
            value: 1000,
        });
    });

    test('reverts required text fields when blurred empty', () => {
        expect(resolveTextFieldBlur('   ', 'Timing A', false)).toEqual({
            action: 'revert',
            error: 'Value required.',
            nextDraft: 'Timing A',
        });
    });

    test('commits non-empty text drafts', () => {
        expect(
            resolveTextFieldBlur('Updated Timing', 'Timing A', false)
        ).toEqual({
            action: 'commit',
            nextDraft: 'Updated Timing',
            value: 'Updated Timing',
        });
    });

    test('formats null numeric values as an empty draft', () => {
        expect(formatNumberFieldValue(null)).toBe('');
    });
});
