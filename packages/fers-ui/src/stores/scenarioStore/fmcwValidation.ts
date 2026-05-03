// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import {
    getDechirpMode,
    isDechirpReferenceSource,
    isFmcwWaveformType,
    isRecord,
} from './fmcwModeConfig';
import type {
    GlobalParameters,
    PlatformComponent,
    ScenarioData,
    SchedulePeriod,
    Waveform,
} from './types';

export type FmcwValidationSeverity = 'error' | 'warning';

export type FmcwValidationIssue = {
    severity: FmcwValidationSeverity;
    message: string;
    itemId?: string;
    componentId?: string;
    waveformId?: string;
    field?: string;
};

type FmcwWaveform = Extract<
    Waveform,
    { waveformType: 'fmcw_linear_chirp' | 'fmcw_triangle' }
>;

type FmcwEmitterComponent = Extract<
    PlatformComponent,
    { type: 'transmitter' | 'monostatic' }
>;

const TRIANGLE_EPSILON = 1e-12;
const IF_CHAIN_FIELD_KEYS = [
    'if_sample_rate',
    'if_filter_bandwidth',
    'if_filter_transition_width',
] as const;

const isFmcwWaveform = (
    waveform: Waveform | undefined
): waveform is FmcwWaveform =>
    waveform?.waveformType === 'fmcw_linear_chirp' ||
    waveform?.waveformType === 'fmcw_triangle';

const formatNumber = (value: number): string =>
    value.toLocaleString(undefined, { maximumSignificantDigits: 6 });

function pushIssue(
    issues: FmcwValidationIssue[],
    issue: FmcwValidationIssue
): void {
    issues.push(issue);
}

function hasIfChainFields(config: unknown): boolean {
    return (
        isRecord(config) &&
        IF_CHAIN_FIELD_KEYS.some((key) => Object.hasOwn(config, key))
    );
}

function getValidIfChainNumber(
    config: unknown,
    key: (typeof IF_CHAIN_FIELD_KEYS)[number],
    label: string,
    component: Pick<PlatformComponent, 'id' | 'name'>,
    issues: FmcwValidationIssue[]
): number | undefined {
    if (!isRecord(config) || !Object.hasOwn(config, key)) {
        return undefined;
    }

    const value = config[key];
    if (typeof value !== 'number' || !Number.isFinite(value) || value <= 0) {
        pushIssue(issues, {
            severity: 'error',
            itemId: component.id,
            componentId: component.id,
            field: 'fmcwModeConfig',
            message: `${component.name} ${label} must be a finite positive value.`,
        });
        return undefined;
    }
    return value;
}

function effectiveSchedule(
    schedule: SchedulePeriod[],
    globalParameters: GlobalParameters
): SchedulePeriod[] {
    return schedule.length > 0
        ? schedule
        : [{ start: globalParameters.start, end: globalParameters.end }];
}

export function validateFmcwWaveform(
    waveform: Waveform,
    globalParameters: GlobalParameters
): FmcwValidationIssue[] {
    if (!isFmcwWaveform(waveform)) {
        return [];
    }

    const issues: FmcwValidationIssue[] = [];
    const sweepStart = waveform.start_frequency_offset ?? 0;
    const sweepEnd =
        waveform.waveformType === 'fmcw_linear_chirp' &&
        waveform.direction === 'down'
            ? sweepStart - waveform.chirp_bandwidth
            : sweepStart + waveform.chirp_bandwidth;
    const fLow = Math.min(sweepStart, sweepEnd);
    const fHigh = Math.max(sweepStart, sweepEnd);
    const maxBaseband = Math.max(Math.abs(fLow), Math.abs(fHigh));
    const effectiveRate =
        globalParameters.rate * globalParameters.oversample_ratio;

    if (
        waveform.waveformType === 'fmcw_linear_chirp' &&
        waveform.chirp_period < waveform.chirp_duration
    ) {
        pushIssue(issues, {
            severity: 'error',
            itemId: waveform.id,
            waveformId: waveform.id,
            field: 'chirp_period',
            message:
                'Chirp period must be greater than or equal to chirp duration.',
        });
    }

    if (effectiveRate <= maxBaseband) {
        pushIssue(issues, {
            severity: 'error',
            itemId: waveform.id,
            waveformId: waveform.id,
            message: `Effective sample rate ${formatNumber(
                effectiveRate
            )} Hz must exceed FMCW sweep baseband ${formatNumber(
                maxBaseband
            )} Hz.`,
        });
    } else if (maxBaseband > 0 && effectiveRate < 1.1 * maxBaseband) {
        pushIssue(issues, {
            severity: 'warning',
            itemId: waveform.id,
            waveformId: waveform.id,
            message: `Effective sample rate ${formatNumber(
                effectiveRate
            )} Hz is within 10% of the FMCW aliasing limit ${formatNumber(
                maxBaseband
            )} Hz.`,
        });
    }

    if (waveform.carrier_frequency + fLow <= 0) {
        pushIssue(issues, {
            severity: 'error',
            itemId: waveform.id,
            waveformId: waveform.id,
            message:
                'Carrier frequency plus the lower sweep edge must stay positive.',
        });
    }

    return issues;
}

function validateFmcwEmitterSchedule(
    component: FmcwEmitterComponent,
    waveform: FmcwWaveform,
    globalParameters: GlobalParameters
): FmcwValidationIssue[] {
    const issues: FmcwValidationIssue[] = [];
    const schedule = effectiveSchedule(component.schedule, globalParameters);

    for (const period of schedule) {
        const duration = period.end - period.start;
        if (waveform.waveformType === 'fmcw_linear_chirp') {
            if (duration < waveform.chirp_duration) {
                pushIssue(issues, {
                    severity: 'error',
                    itemId: component.id,
                    componentId: component.id,
                    waveformId: waveform.id,
                    field: 'schedule',
                    message: `${component.name} has schedule duration ${formatNumber(
                        duration
                    )} s shorter than FMCW chirp duration ${formatNumber(
                        waveform.chirp_duration
                    )} s.`,
                });
            } else if (duration < waveform.chirp_period) {
                pushIssue(issues, {
                    severity: 'warning',
                    itemId: component.id,
                    componentId: component.id,
                    waveformId: waveform.id,
                    field: 'schedule',
                    message: `${component.name} has schedule duration ${formatNumber(
                        duration
                    )} s shorter than FMCW chirp period ${formatNumber(
                        waveform.chirp_period
                    )} s.`,
                });
            }
            continue;
        }

        const trianglePeriod = 2 * waveform.chirp_duration;
        if (duration < trianglePeriod) {
            pushIssue(issues, {
                severity: 'error',
                itemId: component.id,
                componentId: component.id,
                waveformId: waveform.id,
                field: 'schedule',
                message: `${component.name} has schedule duration ${formatNumber(
                    duration
                )} s shorter than FMCW triangle period ${formatNumber(
                    trianglePeriod
                )} s.`,
            });
            continue;
        }

        const fullTriangles = Math.floor(duration / trianglePeriod);
        const leftover = duration - fullTriangles * trianglePeriod;
        if (leftover > TRIANGLE_EPSILON) {
            pushIssue(issues, {
                severity: 'warning',
                itemId: component.id,
                componentId: component.id,
                waveformId: waveform.id,
                field: 'schedule',
                message: `${component.name} schedule leaves ${formatNumber(
                    leftover
                )} s silent after the last complete FMCW triangle.`,
            });
        }
    }

    return issues;
}

function validateFmcwReceiverDechirpConfig(
    component: Extract<PlatformComponent, { type: 'monostatic' | 'receiver' }>,
    fmcwEmitterNames: ReadonlySet<string>,
    fmcwWaveformNames: ReadonlySet<string>
): FmcwValidationIssue[] {
    const issues: FmcwValidationIssue[] = [];

    if (component.radarType !== 'fmcw') {
        return issues;
    }

    const config = component.fmcwModeConfig;
    const mode = getDechirpMode(config);
    const reference =
        isRecord(config) && isRecord(config.dechirp_reference)
            ? config.dechirp_reference
            : null;

    if (mode === 'none') {
        if (reference) {
            pushIssue(issues, {
                severity: 'error',
                itemId: component.id,
                componentId: component.id,
                field: 'fmcwModeConfig',
                message: `${component.name} declares a dechirp reference while dechirp mode is none.`,
            });
        }
        if (hasIfChainFields(config)) {
            pushIssue(issues, {
                severity: 'error',
                itemId: component.id,
                componentId: component.id,
                field: 'fmcwModeConfig',
                message: `${component.name} declares IF-chain settings while dechirp mode is none.`,
            });
        }
        return issues;
    }

    const ifSampleRate = getValidIfChainNumber(
        config,
        'if_sample_rate',
        'IF sample rate',
        component,
        issues
    );
    const ifFilterBandwidth = getValidIfChainNumber(
        config,
        'if_filter_bandwidth',
        'IF filter bandwidth',
        component,
        issues
    );
    getValidIfChainNumber(
        config,
        'if_filter_transition_width',
        'IF transition width',
        component,
        issues
    );
    if (
        ifSampleRate === undefined &&
        (ifFilterBandwidth !== undefined ||
            (isRecord(config) &&
                Object.hasOwn(config, 'if_filter_transition_width')))
    ) {
        pushIssue(issues, {
            severity: 'error',
            itemId: component.id,
            componentId: component.id,
            field: 'fmcwModeConfig',
            message: `${component.name} IF filter settings require an IF sample rate.`,
        });
    }
    if (
        ifSampleRate !== undefined &&
        ifFilterBandwidth !== undefined &&
        ifFilterBandwidth >= ifSampleRate / 2
    ) {
        pushIssue(issues, {
            severity: 'error',
            itemId: component.id,
            componentId: component.id,
            field: 'fmcwModeConfig',
            message: `${component.name} IF filter bandwidth must be less than half the IF sample rate.`,
        });
    }

    if (!reference) {
        pushIssue(issues, {
            severity: 'error',
            itemId: component.id,
            componentId: component.id,
            field: 'fmcwModeConfig',
            message: `${component.name} enables ${mode} dechirping but does not declare a dechirp reference.`,
        });
        return issues;
    }

    if (!isDechirpReferenceSource(reference.source)) {
        pushIssue(issues, {
            severity: 'error',
            itemId: component.id,
            componentId: component.id,
            field: 'fmcwModeConfig',
            message: `${component.name} dechirp reference source must be attached, transmitter, or custom.`,
        });
        return issues;
    }

    switch (reference.source) {
        case 'attached':
            if (component.type !== 'monostatic') {
                pushIssue(issues, {
                    severity: 'error',
                    itemId: component.id,
                    componentId: component.id,
                    field: 'fmcwModeConfig',
                    message: `${component.name} uses an attached dechirp reference, but only monostatic receivers have an attached transmitter.`,
                });
            }
            if (
                'transmitter_name' in reference ||
                'waveform_name' in reference
            ) {
                pushIssue(issues, {
                    severity: 'error',
                    itemId: component.id,
                    componentId: component.id,
                    field: 'fmcwModeConfig',
                    message: `${component.name} attached dechirp reference must not set transmitter or waveform names.`,
                });
            }
            break;
        case 'transmitter': {
            const transmitterName =
                typeof reference.transmitter_name === 'string'
                    ? reference.transmitter_name
                    : '';
            if (transmitterName.trim().length === 0) {
                pushIssue(issues, {
                    severity: 'error',
                    itemId: component.id,
                    componentId: component.id,
                    field: 'fmcwModeConfig',
                    message: `${component.name} transmitter dechirp reference requires a transmitter name.`,
                });
                break;
            }
            if ('waveform_name' in reference) {
                pushIssue(issues, {
                    severity: 'error',
                    itemId: component.id,
                    componentId: component.id,
                    field: 'fmcwModeConfig',
                    message: `${component.name} transmitter dechirp reference must not set a waveform name.`,
                });
            }
            if (!fmcwEmitterNames.has(transmitterName)) {
                pushIssue(issues, {
                    severity: 'error',
                    itemId: component.id,
                    componentId: component.id,
                    field: 'fmcwModeConfig',
                    message: `${component.name} dechirp reference transmitter '${transmitterName}' must be an FMCW transmitter with an FMCW waveform.`,
                });
            }
            break;
        }
        case 'custom': {
            const waveformName =
                typeof reference.waveform_name === 'string'
                    ? reference.waveform_name
                    : '';
            if (waveformName.trim().length === 0) {
                pushIssue(issues, {
                    severity: 'error',
                    itemId: component.id,
                    componentId: component.id,
                    field: 'fmcwModeConfig',
                    message: `${component.name} custom dechirp reference requires a waveform name.`,
                });
                break;
            }
            if ('transmitter_name' in reference) {
                pushIssue(issues, {
                    severity: 'error',
                    itemId: component.id,
                    componentId: component.id,
                    field: 'fmcwModeConfig',
                    message: `${component.name} custom dechirp reference must not set a transmitter name.`,
                });
            }
            if (!fmcwWaveformNames.has(waveformName)) {
                pushIssue(issues, {
                    severity: 'error',
                    itemId: component.id,
                    componentId: component.id,
                    field: 'fmcwModeConfig',
                    message: `${component.name} custom dechirp reference waveform '${waveformName}' must be a top-level FMCW waveform.`,
                });
            }
            break;
        }
    }

    return issues;
}

export function validateFmcwScenario(
    scenario: Pick<ScenarioData, 'globalParameters' | 'waveforms' | 'platforms'>
): FmcwValidationIssue[] {
    const issues = scenario.waveforms.flatMap((waveform) =>
        validateFmcwWaveform(waveform, scenario.globalParameters)
    );
    const waveformsById = new Map(
        scenario.waveforms.map((waveform) => [waveform.id, waveform])
    );
    const fmcwWaveformNames = new Set(
        scenario.waveforms
            .filter((waveform) => isFmcwWaveformType(waveform.waveformType))
            .map((waveform) => waveform.name)
    );
    const fmcwEmitterNames = new Set(
        scenario.platforms.flatMap((platform) =>
            platform.components.flatMap((component) => {
                if (
                    component.type !== 'transmitter' &&
                    component.type !== 'monostatic'
                ) {
                    return [];
                }
                const waveform = component.waveformId
                    ? waveformsById.get(component.waveformId)
                    : undefined;
                return component.radarType === 'fmcw' &&
                    isFmcwWaveform(waveform)
                    ? [component.name]
                    : [];
            })
        )
    );

    for (const platform of scenario.platforms) {
        for (const component of platform.components) {
            if (
                component.type === 'receiver' ||
                component.type === 'monostatic'
            ) {
                issues.push(
                    ...validateFmcwReceiverDechirpConfig(
                        component,
                        fmcwEmitterNames,
                        fmcwWaveformNames
                    )
                );
            }

            if (
                component.type !== 'transmitter' &&
                component.type !== 'monostatic'
            ) {
                continue;
            }

            if (component.radarType !== 'fmcw' || !component.waveformId) {
                continue;
            }

            const waveform = waveformsById.get(component.waveformId);
            if (!isFmcwWaveform(waveform)) {
                pushIssue(issues, {
                    severity: 'error',
                    itemId: component.id,
                    componentId: component.id,
                    waveformId: component.waveformId,
                    message: `${component.name} is FMCW but does not reference an FMCW waveform.`,
                });
                continue;
            }

            issues.push(
                ...validateFmcwEmitterSchedule(
                    component,
                    waveform,
                    scenario.globalParameters
                )
            );
        }
    }

    return issues;
}

export function getBlockingFmcwValidationMessage(
    scenario: Pick<ScenarioData, 'globalParameters' | 'waveforms' | 'platforms'>
): string | null {
    const firstError = validateFmcwScenario(scenario).find(
        (issue) => issue.severity === 'error'
    );
    return firstError?.message ?? null;
}
