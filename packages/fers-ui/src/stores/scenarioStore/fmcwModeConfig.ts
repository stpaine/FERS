// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

export type DechirpMode = 'none' | 'physical' | 'ideal';
export type DechirpReferenceSource = 'attached' | 'transmitter' | 'custom';

export type DechirpReference =
    | { source: 'attached' }
    | { source: 'transmitter'; transmitter_name?: string }
    | { source: 'custom'; waveform_name?: string };

export type ReceiverFmcwModeConfig = {
    dechirp_mode?: DechirpMode;
    dechirp_reference?: DechirpReference;
    if_sample_rate?: number;
    if_filter_bandwidth?: number;
    if_filter_transition_width?: number;
};

export type FmcwIfChainConfig = Pick<
    ReceiverFmcwModeConfig,
    'if_sample_rate' | 'if_filter_bandwidth' | 'if_filter_transition_width'
>;

export const DECHIRP_MODE_OPTIONS: ReadonlyArray<{
    value: DechirpMode;
    label: string;
}> = [
    { value: 'none', label: 'None' },
    { value: 'physical', label: 'Physical' },
    { value: 'ideal', label: 'Ideal' },
];

export const DECHIRP_REFERENCE_SOURCE_OPTIONS: ReadonlyArray<{
    value: DechirpReferenceSource;
    label: string;
}> = [
    { value: 'attached', label: 'Attached' },
    { value: 'transmitter', label: 'Transmitter' },
    { value: 'custom', label: 'Custom Waveform' },
];

const DECHIRP_MODES = new Set<DechirpMode>(['none', 'physical', 'ideal']);
const DECHIRP_REFERENCE_SOURCES = new Set<DechirpReferenceSource>([
    'attached',
    'transmitter',
    'custom',
]);

export const FMCW_WAVEFORM_TYPES = [
    'fmcw_linear_chirp',
    'fmcw_triangle',
] as const;

export type FmcwWaveformType = (typeof FMCW_WAVEFORM_TYPES)[number];

export function isFmcwWaveformType(
    waveformType: string | undefined
): waveformType is FmcwWaveformType {
    return FMCW_WAVEFORM_TYPES.includes(waveformType as FmcwWaveformType);
}

export function isRecord(value: unknown): value is Record<string, unknown> {
    return typeof value === 'object' && value !== null && !Array.isArray(value);
}

export function isDechirpMode(value: unknown): value is DechirpMode {
    return typeof value === 'string' && DECHIRP_MODES.has(value as DechirpMode);
}

export function isDechirpReferenceSource(
    value: unknown
): value is DechirpReferenceSource {
    return (
        typeof value === 'string' &&
        DECHIRP_REFERENCE_SOURCES.has(value as DechirpReferenceSource)
    );
}

export function getDechirpMode(config: unknown): DechirpMode {
    if (!isRecord(config)) {
        return 'none';
    }
    return isDechirpMode(config.dechirp_mode) ? config.dechirp_mode : 'none';
}

export function getDechirpReferenceSource(
    reference: unknown
): DechirpReferenceSource {
    if (!isRecord(reference)) {
        return 'attached';
    }
    return isDechirpReferenceSource(reference.source)
        ? reference.source
        : 'attached';
}

export function getFmcwIfChainConfig(config: unknown): FmcwIfChainConfig {
    if (!isRecord(config)) {
        return {};
    }

    const ifChainConfig: FmcwIfChainConfig = {};
    if (typeof config.if_sample_rate === 'number') {
        ifChainConfig.if_sample_rate = config.if_sample_rate;
    }
    if (typeof config.if_filter_bandwidth === 'number') {
        ifChainConfig.if_filter_bandwidth = config.if_filter_bandwidth;
    }
    if (typeof config.if_filter_transition_width === 'number') {
        ifChainConfig.if_filter_transition_width =
            config.if_filter_transition_width;
    }
    return ifChainConfig;
}

export function createDechirpReference(
    source: DechirpReferenceSource,
    currentReference?: unknown
): DechirpReference {
    const current = isRecord(currentReference) ? currentReference : {};
    switch (source) {
        case 'attached':
            return { source };
        case 'transmitter':
            return {
                source,
                ...(typeof current.transmitter_name === 'string'
                    ? { transmitter_name: current.transmitter_name }
                    : {}),
            };
        case 'custom':
            return {
                source,
                ...(typeof current.waveform_name === 'string'
                    ? { waveform_name: current.waveform_name }
                    : {}),
            };
    }
}

export function normalizeFmcwModeConfig(
    config: unknown
): ReceiverFmcwModeConfig {
    if (!isRecord(config)) {
        return {};
    }

    const dechirpMode = getDechirpMode(config);
    if (dechirpMode === 'none') {
        return {};
    }

    const reference = createDechirpReference(
        getDechirpReferenceSource(config.dechirp_reference),
        config.dechirp_reference
    );

    return {
        dechirp_mode: dechirpMode,
        dechirp_reference: reference,
        ...getFmcwIfChainConfig(config),
    };
}

export function createFmcwModeConfig(
    dechirpMode: DechirpMode,
    currentConfig?: unknown
): ReceiverFmcwModeConfig {
    if (dechirpMode === 'none') {
        return {};
    }

    const current = normalizeFmcwModeConfig(currentConfig);
    return {
        dechirp_mode: dechirpMode,
        dechirp_reference:
            current.dechirp_reference ?? createDechirpReference('attached'),
        ...getFmcwIfChainConfig(current),
    };
}
