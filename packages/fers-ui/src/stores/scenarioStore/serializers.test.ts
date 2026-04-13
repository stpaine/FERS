// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import { serializeAntenna } from './serializers';
import type { Antenna } from './types';

describe('serializeAntenna', () => {
    test('uses an isotropic backend placeholder while an H5 antenna has no filename', () => {
        const antenna: Antenna = {
            id: '1',
            type: 'Antenna',
            name: 'Draft H5',
            pattern: 'file',
            filename: '   ',
            efficiency: 0.8,
            meshScale: 1,
            design_frequency: null,
        };

        expect(serializeAntenna(antenna)).toEqual({
            id: '1',
            name: 'Draft H5',
            pattern: 'isotropic',
            efficiency: 0.8,
        });
    });

    test('preserves H5 antenna payload once a filename is present', () => {
        const antenna: Antenna = {
            id: '2',
            type: 'Antenna',
            name: 'Loaded H5',
            pattern: 'file',
            filename: '/tmp/pattern.h5',
            efficiency: 1,
            meshScale: 1,
            design_frequency: null,
        };

        expect(serializeAntenna(antenna)).toEqual({
            id: '2',
            name: 'Loaded H5',
            pattern: 'file',
            filename: '/tmp/pattern.h5',
            efficiency: 1,
        });
    });
});
