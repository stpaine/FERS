// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import DeleteIcon from '@mui/icons-material/Delete';
import {
    Alert,
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
    Platform,
    PlatformComponent,
    ReceiverComponent,
    SchedulePeriod,
    TargetComponent,
    TransmitterComponent,
    useScenarioStore,
    Waveform,
} from '@/stores/scenarioStore';
import {
    createDechirpReference,
    createFmcwModeConfig,
    DECHIRP_MODE_OPTIONS,
    DECHIRP_REFERENCE_SOURCE_OPTIONS,
    DechirpMode,
    DechirpReferenceSource,
    getDechirpMode,
    getDechirpReferenceSource,
    isFmcwWaveformType,
    ReceiverFmcwModeConfig,
} from '@/stores/scenarioStore/fmcwModeConfig';
import { validateFmcwScenario } from '@/stores/scenarioStore/fmcwValidation';
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

const WAVEFORM_TYPE_BY_RADAR_TYPE: Record<RadarType, string[]> = {
    pulsed: ['pulsed_from_file'],
    cw: ['cw'],
    fmcw: ['fmcw_linear_chirp', 'fmcw_triangle'],
};

export function isWaveformCompatibleWithRadarType(
    waveform: CompatibleWaveform | undefined,
    radarType: RadarType
): boolean {
    return waveform
        ? WAVEFORM_TYPE_BY_RADAR_TYPE[radarType].includes(waveform.waveformType)
        : false;
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

export {
    createDechirpReference,
    createFmcwModeConfig,
    DECHIRP_MODE_OPTIONS,
    DECHIRP_REFERENCE_SOURCE_OPTIONS,
};

const uniqueNames = (names: string[]): string[] =>
    Array.from(new Set(names.filter((name) => name.trim().length > 0)));

const includeCurrentName = (names: string[], currentName?: string): string[] =>
    currentName && !names.includes(currentName)
        ? [currentName, ...names]
        : names;

export function getFmcwWaveformNames(waveforms: Waveform[]): string[] {
    return uniqueNames(
        waveforms
            .filter((waveform) => isFmcwWaveformType(waveform.waveformType))
            .map((waveform) => waveform.name)
    );
}

export function getFmcwEmitterNames(
    platforms: Platform[],
    waveforms: Waveform[]
): string[] {
    const waveformsById = new Map(
        waveforms.map((waveform) => [waveform.id, waveform])
    );
    return uniqueNames(
        platforms.flatMap((platform) =>
            platform.components.flatMap((component) => {
                if (
                    component.type !== 'transmitter' &&
                    component.type !== 'monostatic'
                ) {
                    return [];
                }
                if (
                    component.radarType !== 'fmcw' ||
                    !component.waveformId ||
                    !isFmcwWaveformType(
                        waveformsById.get(component.waveformId)?.waveformType
                    )
                ) {
                    return [];
                }
                return [component.name];
            })
        )
    );
}

export function getAvailableDechirpReferenceSourceOptions(
    componentType: MonostaticComponent['type'] | ReceiverComponent['type'],
    currentSource?: DechirpReferenceSource
) {
    return DECHIRP_REFERENCE_SOURCE_OPTIONS.filter(
        (option) =>
            option.value !== 'attached' ||
            componentType === 'monostatic' ||
            currentSource === 'attached'
    );
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
    const {
        updateItem,
        waveforms,
        timings,
        antennas,
        platforms,
        globalParameters,
        setPlatformRcsModel,
    } = useScenarioStore.getState();
    const fmcwIssues = validateFmcwScenario({
        globalParameters,
        waveforms,
        platforms,
    });
    const fmcwEmitterNames = getFmcwEmitterNames(platforms, waveforms);
    const fmcwWaveformNames = getFmcwWaveformNames(waveforms);

    // Updates are targeted using the array index in the path string
    const handleChange = (path: string, value: unknown) =>
        updateItem(platformId, `components.${index}.${path}`, value);
    const handleComponentChange = (value: PlatformComponent) =>
        updateItem(platformId, `components.${index}`, value);

    const renderSchedule = (
        c: MonostaticComponent | TransmitterComponent | ReceiverComponent
    ) => {
        const schedule = c.schedule || [];
        const scheduleIssues = fmcwIssues.filter(
            (issue) => issue.componentId === c.id && issue.field === 'schedule'
        );

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
                {scheduleIssues.map((issue) => (
                    <Alert
                        key={issue.message}
                        severity={issue.severity}
                        variant="outlined"
                    >
                        {issue.message}
                    </Alert>
                ))}
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
            const prepareModeChange = <
                T extends
                    | MonostaticComponent
                    | TransmitterComponent
                    | ReceiverComponent,
            >(
                component: T,
                waveformId?: string | null
            ): T => {
                const nextComponent = {
                    ...component,
                    radarType: nextRadarType,
                    ...(waveformId !== undefined ? { waveformId } : {}),
                } as T;

                if (
                    nextComponent.type === 'receiver' ||
                    nextComponent.type === 'monostatic'
                ) {
                    if (nextRadarType === 'fmcw') {
                        return {
                            ...nextComponent,
                            fmcwModeConfig: nextComponent.fmcwModeConfig ?? {},
                        } as T;
                    }

                    const {
                        fmcwModeConfig: _fmcwModeConfig,
                        ...withoutFmcwMode
                    } = nextComponent;
                    return withoutFmcwMode as T;
                }

                return nextComponent;
            };

            if ('waveformId' in c) {
                const waveformId = shouldClearWaveformForRadarType(
                    c.waveformId,
                    waveforms,
                    nextRadarType
                )
                    ? null
                    : c.waveformId;

                handleComponentChange(prepareModeChange(c, waveformId));
                return;
            }

            handleComponentChange(prepareModeChange(c));
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
                            handleChange(
                                'antennaId',
                                e.target.value === '' ? null : e.target.value
                            )
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
                            handleChange(
                                'timingId',
                                e.target.value === '' ? null : e.target.value
                            )
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

    const renderFmcwReceiverFields = (
        c: MonostaticComponent | ReceiverComponent
    ) => {
        const config = c.fmcwModeConfig ?? {};
        const dechirpMode = getDechirpMode(config);
        const reference = config.dechirp_reference;
        const referenceSource = getDechirpReferenceSource(reference);
        const transmitterName =
            reference?.source === 'transmitter'
                ? reference.transmitter_name
                : undefined;
        const waveformName =
            reference?.source === 'custom'
                ? reference.waveform_name
                : undefined;
        const referenceIssues = fmcwIssues.filter(
            (issue) =>
                issue.componentId === c.id && issue.field === 'fmcwModeConfig'
        );
        const referenceSourceOptions =
            getAvailableDechirpReferenceSourceOptions(c.type, referenceSource);

        const createDefaultReference = () => {
            if (c.type === 'monostatic') {
                return createDechirpReference('attached');
            }
            const transmitterName = fmcwEmitterNames[0];
            if (transmitterName) {
                return {
                    source: 'transmitter' as const,
                    transmitter_name: transmitterName,
                };
            }
            const waveformName = fmcwWaveformNames[0];
            if (waveformName) {
                return {
                    source: 'custom' as const,
                    waveform_name: waveformName,
                };
            }
            return createDechirpReference('transmitter');
        };

        const commitConfig = (nextConfig: ReceiverFmcwModeConfig) => {
            handleChange('fmcwModeConfig', nextConfig);
        };

        const handleModeChange = (nextMode: DechirpMode) => {
            const nextConfig = createFmcwModeConfig(nextMode, config);
            if (
                nextMode !== 'none' &&
                (!nextConfig.dechirp_reference ||
                    (c.type === 'receiver' &&
                        nextConfig.dechirp_reference.source === 'attached'))
            ) {
                nextConfig.dechirp_reference = createDefaultReference();
            }
            commitConfig(nextConfig);
        };

        const handleReferenceSourceChange = (
            nextSource: DechirpReferenceSource
        ) => {
            const nextReference = createDechirpReference(nextSource, reference);
            if (
                nextReference.source === 'transmitter' &&
                !nextReference.transmitter_name &&
                fmcwEmitterNames[0]
            ) {
                nextReference.transmitter_name = fmcwEmitterNames[0];
            }
            if (
                nextReference.source === 'custom' &&
                !nextReference.waveform_name &&
                fmcwWaveformNames[0]
            ) {
                nextReference.waveform_name = fmcwWaveformNames[0];
            }
            commitConfig({
                dechirp_mode: dechirpMode,
                dechirp_reference: nextReference,
            });
        };

        const handleTransmitterReferenceChange = (name: string) => {
            commitConfig({
                dechirp_mode: dechirpMode,
                dechirp_reference: {
                    source: 'transmitter',
                    ...(name ? { transmitter_name: name } : {}),
                },
            });
        };

        const handleCustomWaveformReferenceChange = (name: string) => {
            commitConfig({
                dechirp_mode: dechirpMode,
                dechirp_reference: {
                    source: 'custom',
                    ...(name ? { waveform_name: name } : {}),
                },
            });
        };

        return (
            <Section title="FMCW Receiver">
                {referenceIssues.map((issue) => (
                    <Alert
                        key={issue.message}
                        severity={issue.severity}
                        variant="outlined"
                    >
                        {issue.message}
                    </Alert>
                ))}
                <FormControl fullWidth size="small">
                    <InputLabel>Dechirp Mode</InputLabel>
                    <Select
                        label="Dechirp Mode"
                        value={dechirpMode}
                        onChange={(e) =>
                            handleModeChange(e.target.value as DechirpMode)
                        }
                    >
                        {DECHIRP_MODE_OPTIONS.map((option) => (
                            <MenuItem key={option.value} value={option.value}>
                                {option.label}
                            </MenuItem>
                        ))}
                    </Select>
                </FormControl>

                {dechirpMode !== 'none' && (
                    <>
                        <FormControl fullWidth size="small">
                            <InputLabel>Dechirp Reference</InputLabel>
                            <Select
                                label="Dechirp Reference"
                                value={referenceSource}
                                onChange={(e) =>
                                    handleReferenceSourceChange(
                                        e.target.value as DechirpReferenceSource
                                    )
                                }
                            >
                                {referenceSourceOptions.map((option) => (
                                    <MenuItem
                                        key={option.value}
                                        value={option.value}
                                    >
                                        {option.label}
                                    </MenuItem>
                                ))}
                            </Select>
                        </FormControl>

                        {referenceSource === 'transmitter' && (
                            <FormControl fullWidth size="small">
                                <InputLabel>Reference Transmitter</InputLabel>
                                <Select
                                    label="Reference Transmitter"
                                    value={transmitterName ?? ''}
                                    onChange={(e) =>
                                        handleTransmitterReferenceChange(
                                            e.target.value
                                        )
                                    }
                                >
                                    <MenuItem value="">
                                        <em>None</em>
                                    </MenuItem>
                                    {includeCurrentName(
                                        fmcwEmitterNames,
                                        transmitterName
                                    ).map((name) => (
                                        <MenuItem key={name} value={name}>
                                            {name}
                                        </MenuItem>
                                    ))}
                                </Select>
                            </FormControl>
                        )}

                        {referenceSource === 'custom' && (
                            <FormControl fullWidth size="small">
                                <InputLabel>Reference Waveform</InputLabel>
                                <Select
                                    label="Reference Waveform"
                                    value={waveformName ?? ''}
                                    onChange={(e) =>
                                        handleCustomWaveformReferenceChange(
                                            e.target.value
                                        )
                                    }
                                >
                                    <MenuItem value="">
                                        <em>None</em>
                                    </MenuItem>
                                    {includeCurrentName(
                                        fmcwWaveformNames,
                                        waveformName
                                    ).map((name) => (
                                        <MenuItem key={name} value={name}>
                                            {name}
                                        </MenuItem>
                                    ))}
                                </Select>
                            </FormControl>
                        )}
                    </>
                )}
            </Section>
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
            {(c.radarType as RadarType) === 'fmcw' &&
                renderFmcwReceiverFields(c)}
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
