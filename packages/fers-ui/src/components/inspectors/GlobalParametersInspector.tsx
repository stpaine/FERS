// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import {
    Box,
    FormControl,
    InputLabel,
    MenuItem,
    Select,
    TextField,
} from '@mui/material';
import { useScenarioStore, GlobalParameters } from '@/stores/scenarioStore';
import { NumberField, Section } from './InspectorControls';

interface GlobalParametersInspectorProps {
    item: GlobalParameters;
}

export function GlobalParametersInspector({
    item,
}: GlobalParametersInspectorProps) {
    const { updateItem } = useScenarioStore.getState();
    const handleChange = (path: string, value: unknown) =>
        updateItem(item.id, path, value);

    return (
        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
            <TextField
                label="Simulation Name"
                variant="outlined"
                size="small"
                fullWidth
                value={item.simulation_name}
                onChange={(e) =>
                    handleChange('simulation_name', e.target.value)
                }
            />
            <NumberField
                label="Start Time (s)"
                value={item.start}
                onChange={(v) => handleChange('start', v)}
            />
            <NumberField
                label="End Time (s)"
                value={item.end}
                onChange={(v) => handleChange('end', v)}
            />
            <NumberField
                label="Output Sampling Rate (Hz)"
                value={item.rate}
                onChange={(v) => handleChange('rate', v)}
            />
            <NumberField
                label="Internal Sim Sampling Rate (Hz)"
                value={item.simSamplingRate}
                onChange={(v) => handleChange('simSamplingRate', v)}
            />
            <NumberField
                label="Speed of Light (m/s)"
                value={item.c}
                onChange={(v) => handleChange('c', v)}
            />
            <NumberField
                label="Random Seed"
                value={item.random_seed}
                onChange={(v) => handleChange('random_seed', v)}
            />
            <NumberField
                label="ADC Bits"
                value={item.adc_bits}
                onChange={(v) => handleChange('adc_bits', v)}
            />
            <NumberField
                label="Oversample Ratio"
                value={item.oversample_ratio}
                onChange={(v) => handleChange('oversample_ratio', v)}
            />

            <Section title="Georeference">
                <NumberField
                    label="Origin Latitude (deg)"
                    value={item.origin.latitude}
                    onChange={(v) => handleChange('origin.latitude', v)}
                />
                <NumberField
                    label="Origin Longitude (deg)"
                    value={item.origin.longitude}
                    onChange={(v) => handleChange('origin.longitude', v)}
                />
                <NumberField
                    label="Origin Altitude (m)"
                    value={item.origin.altitude}
                    onChange={(v) => handleChange('origin.altitude', v)}
                />
                <FormControl fullWidth size="small">
                    <InputLabel>Coordinate System</InputLabel>
                    <Select
                        label="Coordinate System"
                        value={item.coordinateSystem.frame}
                        onChange={(e) =>
                            handleChange(
                                'coordinateSystem.frame',
                                e.target.value
                            )
                        }
                    >
                        <MenuItem value="ENU">ENU (East-North-Up)</MenuItem>
                        <MenuItem value="UTM">UTM</MenuItem>
                        <MenuItem value="ECEF">ECEF</MenuItem>
                    </Select>
                </FormControl>
                {item.coordinateSystem.frame === 'UTM' && (
                    <>
                        <NumberField
                            label="UTM Zone"
                            value={item.coordinateSystem.zone ?? null}
                            onChange={(v) =>
                                handleChange('coordinateSystem.zone', v)
                            }
                        />
                        <FormControl fullWidth size="small">
                            <InputLabel>UTM Hemisphere</InputLabel>
                            <Select
                                label="UTM Hemisphere"
                                value={item.coordinateSystem.hemisphere ?? 'N'}
                                onChange={(e) =>
                                    handleChange(
                                        'coordinateSystem.hemisphere',
                                        e.target.value
                                    )
                                }
                            >
                                <MenuItem value="N">North</MenuItem>
                                <MenuItem value="S">South</MenuItem>
                            </Select>
                        </FormControl>
                    </>
                )}
            </Section>
        </Box>
    );
}
