// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import DeleteIcon from '@mui/icons-material/Delete';
import {
    Box,
    Button,
    Checkbox,
    FormControl,
    FormControlLabel,
    IconButton,
    InputLabel,
    MenuItem,
    Select,
    Typography,
} from '@mui/material';
import {
    MonostaticComponent,
    PlatformComponent,
    ReceiverComponent,
    SchedulePeriod,
    TargetComponent,
    TransmitterComponent,
    useScenarioStore,
} from '@/stores/scenarioStore';
import {
    BufferedTextField,
    FileInput,
    NumberField,
    Section,
} from './InspectorControls';

export type RadarType = 'pulsed' | 'cw' | 'fmcw';
type CompatibleWaveform = {
    id: string;
    name: string;
    waveformType: string;
};

export const RADAR_MODE_OPTIONS: ReadonlyArray<{
    value: RadarType;
    label: string;
}> = [
    { value: 'pulsed', label: 'Pulsed' },
    { value: 'cw', label: 'CW' },
    { value: 'fmcw', label: 'FMCW' },
];

const WAVEFORM_TYPE_BY_RADAR_TYPE: Record<RadarType, string> = {
    pulsed: 'pulsed_from_file',
    cw: 'cw',
    fmcw: 'fmcw_linear_chirp',
};

export function isWaveformCompatibleWithRadarType(
    waveform: CompatibleWaveform | undefined,
    radarType: RadarType
): boolean {
    return waveform?.waveformType === WAVEFORM_TYPE_BY_RADAR_TYPE[radarType];
}

export function getCompatibleWaveforms(
    waveforms: CompatibleWaveform[],
    radarType: RadarType
): CompatibleWaveform[] {
    return waveforms.filter((waveform) =>
        isWaveformCompatibleWithRadarType(waveform, radarType)
    );
}

export function shouldClearWaveformForRadarType(
    waveformId: string | null | undefined,
    waveforms: CompatibleWaveform[],
    radarType: RadarType
): boolean {
    if (!waveformId) {
        return false;
    }

    return !isWaveformCompatibleWithRadarType(
        waveforms.find((waveform) => waveform.id === waveformId),
        radarType
    );
}

export function resolveWaveformSelectValue(
    waveformId: string | null | undefined,
    waveforms: CompatibleWaveform[],
    radarType: RadarType
): string {
    return shouldClearWaveformForRadarType(waveformId, waveforms, radarType)
        ? ''
        : (waveformId ?? '');
}

export function getPulsedRadarFieldLabels(radarType: RadarType): string[] {
    return radarType === 'pulsed'
        ? ['PRF (Hz)', 'Window Skip (s)', 'Window Length (s)']
        : [];
}

interface PlatformComponentInspectorProps {
    component: PlatformComponent;
    platformId: string;
    index: number;
}

export function PlatformComponentInspector({
    component,
    platformId,
    index,
}: PlatformComponentInspectorProps) {
    const { updateItem, waveforms, timings, antennas, setPlatformRcsModel } =
        useScenarioStore.getState();

    // Updates are targeted using the array index in the path string
    const handleChange = (path: string, value: unknown) =>
        updateItem(platformId, `components.${index}.${path}`, value);
    const handleComponentChange = (value: PlatformComponent) =>
        updateItem(platformId, `components.${index}`, value);

    const renderSchedule = (
        c: MonostaticComponent | TransmitterComponent | ReceiverComponent
    ) => {
        const schedule = c.schedule || [];

        const handleAddPeriod = () => {
            handleChange('schedule', [...schedule, { start: 0, end: 0 }]);
        };

        const handleRemovePeriod = (idx: number) => {
            const newSchedule = [...schedule];
            newSchedule.splice(idx, 1);
            handleChange('schedule', newSchedule);
        };

        const handlePeriodChange = (
            idx: number,
            field: keyof SchedulePeriod,
            val: number | null
        ) => {
            const newSchedule = [...schedule];
            newSchedule[idx] = { ...newSchedule[idx], [field]: val ?? 0 };
            handleChange('schedule', newSchedule);
        };

        return (
            <Section title="Operating Schedule">
                {schedule.length === 0 && (
                    <Typography variant="body2" color="text.secondary">
                        No specific schedule defined (always active).
                    </Typography>
                )}
                {schedule.map((period, i) => (
                    <Box
                        key={i}
                        sx={{
                            display: 'flex',
                            alignItems: 'center',
                            gap: 1,
                            p: 1,
                            border: 1,
                            borderColor: 'divider',
                            borderRadius: 1,
                        }}
                    >
                        <NumberField
                            label="Start (s)"
                            value={period.start}
                            emptyBehavior="revert"
                            onChange={(v) => handlePeriodChange(i, 'start', v)}
                        />
                        <NumberField
                            label="End (s)"
                            value={period.end}
                            emptyBehavior="revert"
                            onChange={(v) => handlePeriodChange(i, 'end', v)}
                        />
                        <IconButton
                            size="small"
                            onClick={() => handleRemovePeriod(i)}
                            color="error"
                        >
                            <DeleteIcon fontSize="small" />
                        </IconButton>
                    </Box>
                ))}
                <Button
                    onClick={handleAddPeriod}
                    size="small"
                    variant="outlined"
                    sx={{ mt: 1 }}
                >
                    Add Schedule Period
                </Button>
            </Section>
        );
    };

    const renderCommonRadarFields = (
        c: MonostaticComponent | TransmitterComponent | ReceiverComponent
    ) => {
        const radarType = c.radarType as RadarType;
        const compatibleWaveforms = getCompatibleWaveforms(
            waveforms,
            radarType
        );
        const handleRadarTypeChange = (nextRadarType: RadarType) => {
            if ('waveformId' in c) {
                const waveformId = shouldClearWaveformForRadarType(
                    c.waveformId,
                    waveforms,
                    nextRadarType
                )
                    ? null
                    : c.waveformId;

                handleComponentChange({
                    ...c,
                    radarType: nextRadarType,
                    waveformId,
                });
                return;
            }

            handleComponentChange({
                ...c,
                radarType: nextRadarType,
            });
        };

        return (
            <>
                <BufferedTextField
                    label="Component Name"
                    size="small"
                    fullWidth
                    value={c.name}
                    allowEmpty={false}
                    onChange={(v) => handleChange('name', v)}
                    sx={{ mb: 2 }}
                />
                <FormControl fullWidth size="small" sx={{ mb: 2 }}>
                    <InputLabel>Radar Mode</InputLabel>
                    <Select
                        label="Radar Mode"
                        value={radarType}
                        onChange={(e) =>
                            handleRadarTypeChange(e.target.value as RadarType)
                        }
                    >
                        {RADAR_MODE_OPTIONS.map((option) => (
                            <MenuItem key={option.value} value={option.value}>
                                {option.label}
                            </MenuItem>
                        ))}
                    </Select>
                </FormControl>

                {'waveformId' in c && (
                    <FormControl fullWidth size="small" sx={{ mb: 2 }}>
                        <InputLabel>Waveform</InputLabel>
                        <Select
                            label="Waveform"
                            value={resolveWaveformSelectValue(
                                c.waveformId,
                                waveforms,
                                radarType
                            )}
                            onChange={(e) =>
                                handleChange(
                                    'waveformId',
                                    e.target.value === ''
                                        ? null
                                        : e.target.value
                                )
                            }
                        >
                            <MenuItem value="">
                                <em>None</em>
                            </MenuItem>
                            {compatibleWaveforms.map((w) => (
                                <MenuItem key={w.id} value={w.id}>
                                    {w.name}
                                </MenuItem>
                            ))}
                        </Select>
                    </FormControl>
                )}

                <FormControl fullWidth size="small" sx={{ mb: 2 }}>
                    <InputLabel>Antenna</InputLabel>
                    <Select
                        label="Antenna"
                        value={c.antennaId ?? ''}
                        onChange={(e) =>
                            handleChange('antennaId', e.target.value)
                        }
                    >
                        <MenuItem value="">
                            <em>None</em>
                        </MenuItem>
                        {antennas.map((a) => (
                            <MenuItem key={a.id} value={a.id}>
                                {a.name}
                            </MenuItem>
                        ))}
                    </Select>
                </FormControl>
                <FormControl fullWidth size="small" sx={{ mb: 2 }}>
                    <InputLabel>Timing Source</InputLabel>
                    <Select
                        label="Timing Source"
                        value={c.timingId ?? ''}
                        onChange={(e) =>
                            handleChange('timingId', e.target.value)
                        }
                    >
                        <MenuItem value="">
                            <em>None</em>
                        </MenuItem>
                        {timings.map((t) => (
                            <MenuItem key={t.id} value={t.id}>
                                {t.name}
                            </MenuItem>
                        ))}
                    </Select>
                </FormControl>
            </>
        );
    };

    const renderReceiverFields = (
        c: MonostaticComponent | ReceiverComponent
    ) => (
        <>
            {(c.radarType as RadarType) === 'pulsed' && (
                <>
                    <NumberField
                        label="Window Skip (s)"
                        value={c.window_skip}
                        emptyBehavior="revert"
                        onChange={(v) => handleChange('window_skip', v)}
                    />
                    <NumberField
                        label="Window Length (s)"
                        value={c.window_length}
                        emptyBehavior="revert"
                        onChange={(v) => handleChange('window_length', v)}
                    />
                </>
            )}
            <NumberField
                label="Noise Temperature (K)"
                value={c.noiseTemperature}
                emptyBehavior="revert"
                onChange={(v) => handleChange('noiseTemperature', v)}
            />
            <FormControlLabel
                control={
                    <Checkbox
                        checked={c.noDirectPaths}
                        onChange={(e) =>
                            handleChange('noDirectPaths', e.target.checked)
                        }
                    />
                }
                label="Ignore Direct Paths"
            />
            <FormControlLabel
                control={
                    <Checkbox
                        checked={c.noPropagationLoss}
                        onChange={(e) =>
                            handleChange('noPropagationLoss', e.target.checked)
                        }
                    />
                }
                label="Ignore Propagation Loss"
            />
        </>
    );

    switch (component.type) {
        case 'monostatic':
            return (
                <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
                    {renderCommonRadarFields(component)}
                    {(component.radarType as RadarType) === 'pulsed' && (
                        <NumberField
                            label="PRF (Hz)"
                            value={component.prf}
                            emptyBehavior="revert"
                            onChange={(v) => handleChange('prf', v)}
                        />
                    )}
                    {renderReceiverFields(component)}
                    {renderSchedule(component)}
                </Box>
            );
        case 'transmitter':
            return (
                <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
                    {renderCommonRadarFields(component)}
                    {(component.radarType as RadarType) === 'pulsed' && (
                        <NumberField
                            label="PRF (Hz)"
                            value={component.prf}
                            emptyBehavior="revert"
                            onChange={(v) => handleChange('prf', v)}
                        />
                    )}
                    {renderSchedule(component)}
                </Box>
            );
        case 'receiver':
            return (
                <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
                    {renderCommonRadarFields(component)}
                    {(component.radarType as RadarType) === 'pulsed' && (
                        <NumberField
                            label="PRF (Hz)"
                            value={component.prf}
                            emptyBehavior="revert"
                            onChange={(v) => handleChange('prf', v)}
                        />
                    )}
                    {renderReceiverFields(component)}
                    {renderSchedule(component)}
                </Box>
            );
        case 'target':
            return (
                <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
                    <BufferedTextField
                        label="Component Name"
                        size="small"
                        fullWidth
                        value={component.name}
                        allowEmpty={false}
                        onChange={(v) => handleChange('name', v)}
                    />
                    <FormControl fullWidth size="small">
                        <InputLabel>RCS Type</InputLabel>
                        <Select
                            label="RCS Type"
                            value={component.rcs_type}
                            onChange={(e) =>
                                handleChange('rcs_type', e.target.value)
                            }
                        >
                            <MenuItem value="isotropic">Isotropic</MenuItem>
                            <MenuItem value="file">File</MenuItem>
                        </Select>
                    </FormControl>
                    {component.rcs_type === 'isotropic' && (
                        <NumberField
                            label="RCS Value (m^2)"
                            value={component.rcs_value ?? 0}
                            emptyBehavior="revert"
                            onChange={(v) => handleChange('rcs_value', v)}
                        />
                    )}
                    {component.rcs_type === 'file' && (
                        <FileInput
                            label="RCS File"
                            value={component.rcs_filename}
                            onChange={(v) => handleChange('rcs_filename', v)}
                            filters={[{ name: 'RCS Data', extensions: ['*'] }]}
                        />
                    )}

                    <FormControl fullWidth size="small">
                        <InputLabel>RCS Model</InputLabel>
                        <Select
                            label="RCS Model"
                            value={component.rcs_model}
                            onChange={(e) =>
                                setPlatformRcsModel(
                                    platformId,
                                    component.id,
                                    e.target
                                        .value as TargetComponent['rcs_model']
                                )
                            }
                        >
                            <MenuItem value="constant">Constant</MenuItem>
                            <MenuItem value="chisquare">Chi-Square</MenuItem>
                            <MenuItem value="gamma">Gamma</MenuItem>
                        </Select>
                    </FormControl>
                    {(component.rcs_model === 'chisquare' ||
                        component.rcs_model === 'gamma') && (
                        <NumberField
                            label="K Value"
                            value={component.rcs_k ?? 0}
                            emptyBehavior="revert"
                            onChange={(v) => handleChange('rcs_k', v)}
                        />
                    )}
                </Box>
            );
        default:
            return (
                <Typography color="text.secondary">
                    Unknown component type.
                </Typography>
            );
    }
}
