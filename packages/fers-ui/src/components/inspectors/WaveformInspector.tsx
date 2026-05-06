// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import {
    Alert,
    Box,
    FormControl,
    InputLabel,
    MenuItem,
    Select,
} from '@mui/material';
import { useScenarioStore, Waveform } from '@/stores/scenarioStore';
import { createWaveformForType as createWaveformDefaultsForType } from '@/stores/scenarioStore/defaults';
import { validateFmcwWaveform } from '@/stores/scenarioStore/fmcwValidation';
import { BufferedTextField, FileInput, NumberField } from './InspectorControls';

export type WaveformType =
    | 'pulsed_from_file'
    | 'cw'
    | 'fmcw_linear_chirp'
    | 'fmcw_triangle';

type FmcwLinearChirpFields = {
    direction: 'up' | 'down';
    chirp_bandwidth: number;
    chirp_duration: number;
    chirp_period: number;
    start_frequency_offset: number | null;
    chirp_count: number | null;
};

type FmcwTriangleFields = {
    chirp_bandwidth: number;
    chirp_duration: number;
    start_frequency_offset: number | null;
    triangle_count: number | null;
};

type AuthorableWaveform = Omit<Waveform, 'waveformType'> &
    Partial<FmcwLinearChirpFields> &
    Partial<FmcwTriangleFields> & {
        waveformType: WaveformType;
        filename?: string;
    };

export const WAVEFORM_TYPE_OPTIONS: ReadonlyArray<{
    value: WaveformType;
    label: string;
}> = [
    { value: 'pulsed_from_file', label: 'Pulse File' },
    { value: 'cw', label: 'CW' },
    { value: 'fmcw_linear_chirp', label: 'FMCW Linear Chirp' },
    { value: 'fmcw_triangle', label: 'FMCW Triangle' },
];

const DEFAULT_FMCW_LINEAR_CHIRP_FIELDS =
    createWaveformDefaultsForType('fmcw_linear_chirp');
const DEFAULT_FMCW_TRIANGLE_FIELDS =
    createWaveformDefaultsForType('fmcw_triangle');

interface WaveformInspectorProps {
    item: Waveform;
}

const asNumberOrDefault = (value: unknown, fallback: number): number =>
    typeof value === 'number' && Number.isFinite(value) ? value : fallback;

const asNullableNumber = (value: unknown): number | null =>
    typeof value === 'number' && Number.isFinite(value) ? value : null;

const asNullableNumberOrDefault = (
    value: unknown,
    fallback: number | null
): number | null =>
    typeof value === 'number' && Number.isFinite(value) ? value : fallback;

const asNullableInteger = (value: unknown): number | null =>
    typeof value === 'number' && Number.isFinite(value)
        ? Math.trunc(value)
        : null;

export function createWaveformForType(
    waveform: AuthorableWaveform,
    waveformType: WaveformType
): AuthorableWaveform {
    const nextWaveform = {
        ...waveform,
        waveformType,
    };

    if (waveformType === 'pulsed_from_file') {
        return {
            ...nextWaveform,
            filename:
                typeof waveform.filename === 'string' ? waveform.filename : '',
        };
    }

    if (waveformType === 'fmcw_linear_chirp') {
        return {
            ...nextWaveform,
            direction: waveform.direction === 'down' ? 'down' : 'up',
            chirp_bandwidth: asNumberOrDefault(
                waveform.chirp_bandwidth,
                DEFAULT_FMCW_LINEAR_CHIRP_FIELDS.chirp_bandwidth
            ),
            chirp_duration: asNumberOrDefault(
                waveform.chirp_duration,
                DEFAULT_FMCW_LINEAR_CHIRP_FIELDS.chirp_duration
            ),
            chirp_period: asNumberOrDefault(
                waveform.chirp_period,
                DEFAULT_FMCW_LINEAR_CHIRP_FIELDS.chirp_period
            ),
            start_frequency_offset: asNullableNumberOrDefault(
                waveform.start_frequency_offset,
                DEFAULT_FMCW_LINEAR_CHIRP_FIELDS.start_frequency_offset
            ),
            chirp_count: asNullableInteger(waveform.chirp_count),
        };
    }

    if (waveformType === 'fmcw_triangle') {
        return {
            ...nextWaveform,
            chirp_bandwidth: asNumberOrDefault(
                waveform.chirp_bandwidth,
                DEFAULT_FMCW_TRIANGLE_FIELDS.chirp_bandwidth
            ),
            chirp_duration: asNumberOrDefault(
                waveform.chirp_duration,
                DEFAULT_FMCW_TRIANGLE_FIELDS.chirp_duration
            ),
            start_frequency_offset: asNullableNumberOrDefault(
                waveform.start_frequency_offset,
                DEFAULT_FMCW_TRIANGLE_FIELDS.start_frequency_offset
            ),
            triangle_count: asNullableInteger(waveform.triangle_count),
        };
    }

    return nextWaveform;
}

export function getVisibleWaveformFieldLabels(
    waveformType: WaveformType
): string[] {
    if (waveformType === 'pulsed_from_file') {
        return ['Waveform File (.csv, .h5)'];
    }

    if (waveformType === 'fmcw_linear_chirp') {
        return [
            'Direction',
            'Chirp Bandwidth (Hz)',
            'Chirp Duration (s)',
            'Chirp Period (s)',
            'Start Frequency Offset (Hz)',
            'Chirp Count',
        ];
    }

    if (waveformType === 'fmcw_triangle') {
        return [
            'Chirp Bandwidth (Hz)',
            'Chirp Duration (s)',
            'Start Frequency Offset (Hz)',
            'Triangle Count',
        ];
    }

    return [];
}

export function WaveformInspector({ item }: WaveformInspectorProps) {
    const { updateItem, globalParameters } = useScenarioStore.getState();
    const waveform = item as AuthorableWaveform;
    const fmcwIssues = validateFmcwWaveform(item, globalParameters);
    const getFieldIssue = (field: string) =>
        fmcwIssues.find((issue) => issue.field === field);
    const globalFmcwIssues = fmcwIssues.filter((issue) => !issue.field);
    const handleChange = (path: string, value: unknown) =>
        updateItem(item.id, path, value);
    const handleTypeChange = (waveformType: WaveformType) => {
        const nextWaveform = createWaveformForType(waveform, waveformType);
        const currentValues = waveform as Record<string, unknown>;
        const nextEntries = Object.entries(nextWaveform).filter(
            ([key]) => key !== 'id' && key !== 'type'
        );

        for (const [key, value] of nextEntries.filter(
            ([key]) => key !== 'waveformType'
        )) {
            if (!Object.is(currentValues[key], value)) {
                handleChange(key, value);
            }
        }

        if (waveform.waveformType !== waveformType) {
            handleChange('waveformType', waveformType);
        }
    };

    return (
        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
            <BufferedTextField
                label="Name"
                variant="outlined"
                size="small"
                fullWidth
                value={item.name}
                allowEmpty={false}
                onChange={(v) => handleChange('name', v)}
            />
            <FormControl fullWidth size="small">
                <InputLabel>Type</InputLabel>
                <Select
                    label="Type"
                    value={waveform.waveformType}
                    onChange={(e) =>
                        handleTypeChange(e.target.value as WaveformType)
                    }
                >
                    {WAVEFORM_TYPE_OPTIONS.map((option) => (
                        <MenuItem key={option.value} value={option.value}>
                            {option.label}
                        </MenuItem>
                    ))}
                </Select>
            </FormControl>
            <NumberField
                label="Power (W)"
                value={item.power}
                emptyBehavior="revert"
                onChange={(v) => handleChange('power', v)}
            />
            <NumberField
                label="Carrier Frequency (Hz)"
                value={item.carrier_frequency}
                emptyBehavior="revert"
                onChange={(v) => handleChange('carrier_frequency', v)}
            />
            {globalFmcwIssues.map((issue) => (
                <Alert
                    key={issue.message}
                    severity={issue.severity}
                    variant="outlined"
                >
                    {issue.message}
                </Alert>
            ))}
            {waveform.waveformType === 'pulsed_from_file' && (
                <FileInput
                    label="Waveform File (.csv, .h5)"
                    value={waveform.filename}
                    onChange={(v) => handleChange('filename', v)}
                    filters={[
                        { name: 'Waveform', extensions: ['csv', 'h5'] },
                        { name: 'All Files', extensions: ['*'] },
                    ]}
                />
            )}
            {waveform.waveformType === 'fmcw_linear_chirp' && (
                <>
                    <FormControl fullWidth size="small">
                        <InputLabel>Direction</InputLabel>
                        <Select
                            label="Direction"
                            value={waveform.direction ?? 'up'}
                            onChange={(e) =>
                                handleChange('direction', e.target.value)
                            }
                        >
                            <MenuItem value="up">Up</MenuItem>
                            <MenuItem value="down">Down</MenuItem>
                        </Select>
                    </FormControl>
                    <NumberField
                        label="Chirp Bandwidth (Hz)"
                        value={asNumberOrDefault(
                            waveform.chirp_bandwidth,
                            DEFAULT_FMCW_LINEAR_CHIRP_FIELDS.chirp_bandwidth
                        )}
                        emptyBehavior="revert"
                        onChange={(v) => handleChange('chirp_bandwidth', v)}
                    />
                    <NumberField
                        label="Chirp Duration (s)"
                        value={asNumberOrDefault(
                            waveform.chirp_duration,
                            DEFAULT_FMCW_LINEAR_CHIRP_FIELDS.chirp_duration
                        )}
                        emptyBehavior="revert"
                        onChange={(v) => handleChange('chirp_duration', v)}
                    />
                    <NumberField
                        label="Chirp Period (s)"
                        value={asNumberOrDefault(
                            waveform.chirp_period,
                            DEFAULT_FMCW_LINEAR_CHIRP_FIELDS.chirp_period
                        )}
                        emptyBehavior="revert"
                        helperText={getFieldIssue('chirp_period')?.message}
                        externalError={
                            getFieldIssue('chirp_period')?.severity === 'error'
                        }
                        onChange={(v) => handleChange('chirp_period', v)}
                    />
                    <NumberField
                        label="Start Frequency Offset (Hz)"
                        value={asNullableNumber(
                            waveform.start_frequency_offset
                        )}
                        emptyBehavior="null"
                        onChange={(v) =>
                            handleChange('start_frequency_offset', v)
                        }
                    />
                    <NumberField
                        label="Chirp Count"
                        value={asNullableInteger(waveform.chirp_count)}
                        emptyBehavior="null"
                        onChange={(v) =>
                            handleChange(
                                'chirp_count',
                                v === null ? null : Math.trunc(v)
                            )
                        }
                    />
                </>
            )}
            {waveform.waveformType === 'fmcw_triangle' && (
                <>
                    <NumberField
                        label="Chirp Bandwidth (Hz)"
                        value={asNumberOrDefault(
                            waveform.chirp_bandwidth,
                            DEFAULT_FMCW_TRIANGLE_FIELDS.chirp_bandwidth
                        )}
                        emptyBehavior="revert"
                        onChange={(v) => handleChange('chirp_bandwidth', v)}
                    />
                    <NumberField
                        label="Chirp Duration (s)"
                        value={asNumberOrDefault(
                            waveform.chirp_duration,
                            DEFAULT_FMCW_TRIANGLE_FIELDS.chirp_duration
                        )}
                        emptyBehavior="revert"
                        onChange={(v) => handleChange('chirp_duration', v)}
                    />
                    <NumberField
                        label="Start Frequency Offset (Hz)"
                        value={asNullableNumber(
                            waveform.start_frequency_offset
                        )}
                        emptyBehavior="null"
                        onChange={(v) =>
                            handleChange('start_frequency_offset', v)
                        }
                    />
                    <NumberField
                        label="Triangle Count"
                        value={asNullableInteger(waveform.triangle_count)}
                        emptyBehavior="null"
                        onChange={(v) =>
                            handleChange(
                                'triangle_count',
                                v === null ? null : Math.trunc(v)
                            )
                        }
                    />
                </>
            )}
        </Box>
    );
}
