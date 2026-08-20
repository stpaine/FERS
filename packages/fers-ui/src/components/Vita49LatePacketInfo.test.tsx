// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { describe, expect, test } from 'bun:test';
import { renderToStaticMarkup } from 'react-dom/server';
import {
    latePacketTooltipText,
    totalLatePacketCount,
    Vita49LatePacketInfo,
} from './Vita49LatePacketInfo';

describe('VITA49 late-packet information', () => {
    test('explains the deadline and classified packet counts', () => {
        expect(totalLatePacketCount(27, 182)).toBe(209);
        expect(latePacketTooltipText(27, 182)).toBe(
            'Packets dispatched more than 1 ms after their scheduled wall-clock deadline. Data: 27. Context: 182.'
        );
    });

    test('exposes the tooltip explanation to keyboard and screen-reader users', () => {
        const html = renderToStaticMarkup(
            <Vita49LatePacketInfo
                dataPacketCount={27}
                contextPacketCount={182}
            />
        );

        expect(html).toContain('aria-label=');
        expect(html).toContain('Data: 27. Context: 182.');
        expect(html).toContain('<button');
    });
});
