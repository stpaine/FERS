// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

export type SimObjectType =
    | 'Platform'
    | 'Transmitter'
    | 'Receiver'
    | 'Target'
    | 'Antenna'
    | 'Waveform'
    | 'Timing';

const TYPE_CODES: Record<SimObjectType, bigint> = {
    Platform: 1n,
    Transmitter: 2n,
    Receiver: 3n,
    Target: 4n,
    Antenna: 5n,
    Waveform: 6n,
    Timing: 7n,
};

const CODE_TO_TYPE: Record<number, SimObjectType> = {
    1: 'Platform',
    2: 'Transmitter',
    3: 'Receiver',
    4: 'Target',
    5: 'Antenna',
    6: 'Waveform',
    7: 'Timing',
};

const MAX_COUNTER = 0x0000ffffffffffffn;

const counters: Record<SimObjectType, bigint> = {
    Platform: 1n,
    Transmitter: 1n,
    Receiver: 1n,
    Target: 1n,
    Antenna: 1n,
    Waveform: 1n,
    Timing: 1n,
};

export const normalizeSimId = (value: unknown): string | null => {
    if (typeof value === 'string') {
        return /^\d+$/.test(value) ? value : null;
    }
    if (typeof value === 'number' && Number.isFinite(value)) {
        return Math.trunc(value).toString();
    }
    return null;
};

export const reserveSimId = (id: string): void => {
    const parsed = BigInt(id);
    const typeCode = Number(parsed >> 48n);
    const counter = parsed & MAX_COUNTER;
    const type = CODE_TO_TYPE[typeCode];
    if (!type) return;
    const next = counter + 1n;
    if (next > counters[type]) {
        counters[type] = next;
    }
};

export const seedSimIdCounters = (
    ids: Array<string | null | undefined>
): void => {
    ids.forEach((id) => {
        if (id) {
            reserveSimId(id);
        }
    });
};

export const generateSimId = (type: SimObjectType): string => {
    const counter = counters[type];
    if (counter > MAX_COUNTER) {
        throw new Error(`SimId counter overflow for ${type}.`);
    }
    const id = (TYPE_CODES[type] << 48n) | counter;
    counters[type] = counter + 1n;
    return id.toString();
};
