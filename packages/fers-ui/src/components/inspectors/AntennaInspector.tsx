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
import { useScenarioStore, Antenna } from '@/stores/scenarioStore';
import { FileInput, NumberField } from './InspectorControls';

interface AntennaInspectorProps {
    item: Antenna;
}

export function AntennaInspector({ item }: AntennaInspectorProps) {
    const { updateItem, setAntennaPattern } = useScenarioStore.getState();
    const handleChange = (path: string, value: unknown) =>
        updateItem(item.id, path, value);

    return (
        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
            <TextField
                label="Name"
                variant="outlined"
                size="small"
                fullWidth
                value={item.name}
                onChange={(e) => handleChange('name', e.target.value)}
            />
            <FormControl fullWidth size="small">
                <InputLabel>Pattern</InputLabel>
                <Select
                    label="Pattern"
                    value={item.pattern}
                    onChange={(e) =>
                        setAntennaPattern(
                            item.id,
                            e.target.value as Antenna['pattern']
                        )
                    }
                >
                    <MenuItem value="isotropic">Isotropic</MenuItem>
                    <MenuItem value="sinc">Sinc</MenuItem>
                    <MenuItem value="gaussian">Gaussian</MenuItem>
                    <MenuItem value="squarehorn">Square Horn</MenuItem>
                    <MenuItem value="parabolic">Parabolic</MenuItem>
                    <MenuItem value="xml">XML</MenuItem>
                    <MenuItem value="file">File (H5)</MenuItem>
                </Select>
            </FormControl>
            <NumberField
                label="Efficiency"
                value={item.efficiency}
                onChange={(v) => handleChange('efficiency', v)}
            />
            <NumberField
                label="Mesh Scale Multiplier"
                value={item.meshScale ?? null}
                onChange={(v) => handleChange('meshScale', v)}
            />

            {item.pattern === 'sinc' && (
                <>
                    <NumberField
                        label="Alpha"
                        value={item.alpha ?? null}
                        onChange={(v) => handleChange('alpha', v)}
                    />
                    <NumberField
                        label="Beta"
                        value={item.beta ?? null}
                        onChange={(v) => handleChange('beta', v)}
                    />
                    <NumberField
                        label="Gamma"
                        value={item.gamma ?? null}
                        onChange={(v) => handleChange('gamma', v)}
                    />
                </>
            )}
            {item.pattern === 'gaussian' && (
                <>
                    <NumberField
                        label="Azimuth Scale"
                        value={item.azscale ?? null}
                        onChange={(v) => handleChange('azscale', v)}
                    />
                    <NumberField
                        label="Elevation Scale"
                        value={item.elscale ?? null}
                        onChange={(v) => handleChange('elscale', v)}
                    />
                </>
            )}
            {(item.pattern === 'squarehorn' ||
                item.pattern === 'parabolic') && (
                <>
                    <NumberField
                        label="Diameter (m)"
                        value={item.diameter ?? null}
                        onChange={(v) => handleChange('diameter', v)}
                    />
                    <NumberField
                        label="Design Frequency (Hz)"
                        value={item.design_frequency ?? null}
                        onChange={(v) => handleChange('design_frequency', v)}
                    />
                </>
            )}

            {(item.pattern === 'xml' || item.pattern === 'file') && (
                <FileInput
                    label="Pattern File"
                    value={item.filename}
                    onChange={(v) => handleChange('filename', v)}
                    filters={[
                        {
                            name: 'Antenna Pattern',
                            extensions:
                                item.pattern === 'xml' ? ['xml'] : ['h5'],
                        },
                        { name: 'All Files', extensions: ['*'] },
                    ]}
                />
            )}
        </Box>
    );
}
