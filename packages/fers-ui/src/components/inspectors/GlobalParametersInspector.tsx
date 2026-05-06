// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import {
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogContentText,
    DialogTitle,
    FormControl,
    InputLabel,
    MenuItem,
    Select,
} from '@mui/material';
import { useState } from 'react';
import { GlobalParameters, useScenarioStore } from '@/stores/scenarioStore';
import { deriveUtmZone } from '@/stores/scenarioStore/serializers';
import { BufferedTextField, NumberField, Section } from './InspectorControls';

interface GlobalParametersInspectorProps {
    item: GlobalParameters;
}

export function GlobalParametersInspector({
    item,
}: GlobalParametersInspectorProps) {
    const updateItem = useScenarioStore((s) => s.updateItem);
    const setRotationAngleUnit = useScenarioStore(
        (s) => s.setRotationAngleUnit
    );
    const platforms = useScenarioStore((s) => s.platforms);
    const handleChange = (path: string, value: unknown) =>
        updateItem(item.id, path, value);
    const [pendingRotationUnit, setPendingRotationUnit] = useState<
        GlobalParameters['rotationAngleUnit'] | null
    >(null);

    const hasExistingRotationValues = platforms.some((platform) => {
        if (platform.rotation.type === 'fixed') {
            return [
                platform.rotation.startAzimuth,
                platform.rotation.startElevation,
                platform.rotation.azimuthRate,
                platform.rotation.elevationRate,
            ].some((value) => value !== 0);
        }
        return platform.rotation.waypoints.some(
            (waypoint) => waypoint.azimuth !== 0 || waypoint.elevation !== 0
        );
    });

    const handleRotationUnitChange = (
        nextUnit: GlobalParameters['rotationAngleUnit']
    ) => {
        if (nextUnit === item.rotationAngleUnit) {
            return;
        }
        if (hasExistingRotationValues) {
            setPendingRotationUnit(nextUnit);
            return;
        }
        setRotationAngleUnit(nextUnit, false);
    };

    const handleCoordinateSystemFrameChange = (
        frame: GlobalParameters['coordinateSystem']['frame']
    ) => {
        if (frame === 'UTM') {
            handleChange('coordinateSystem', {
                frame,
                zone:
                    item.coordinateSystem.zone ??
                    deriveUtmZone(item.origin.longitude),
                hemisphere:
                    item.coordinateSystem.hemisphere ??
                    (item.origin.latitude < 0 ? 'S' : 'N'),
            });
            return;
        }

        handleChange('coordinateSystem', { frame });
    };

    return (
        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
            <FormControl fullWidth size="small">
                <InputLabel>Rotation Angle Unit</InputLabel>
                <Select
                    label="Rotation Angle Unit"
                    value={item.rotationAngleUnit}
                    onChange={(e) =>
                        handleRotationUnitChange(
                            e.target
                                .value as GlobalParameters['rotationAngleUnit']
                        )
                    }
                >
                    <MenuItem value="deg">Degrees</MenuItem>
                    <MenuItem value="rad">Radians</MenuItem>
                </Select>
            </FormControl>
            <BufferedTextField
                label="Simulation Name"
                variant="outlined"
                size="small"
                fullWidth
                value={item.simulation_name}
                allowEmpty={false}
                onChange={(v) => handleChange('simulation_name', v)}
            />
            <NumberField
                label="Start Time (s)"
                value={item.start}
                emptyBehavior="revert"
                onChange={(v) => handleChange('start', v)}
            />
            <NumberField
                label="End Time (s)"
                value={item.end}
                emptyBehavior="revert"
                onChange={(v) => handleChange('end', v)}
            />
            <NumberField
                label="Output Sampling Rate (Hz)"
                value={item.rate}
                emptyBehavior="revert"
                onChange={(v) => handleChange('rate', v)}
            />
            <NumberField
                label="Internal Sim Sampling Rate (Hz)"
                value={item.simSamplingRate}
                emptyBehavior="revert"
                onChange={(v) => handleChange('simSamplingRate', v)}
            />
            <NumberField
                label="Speed of Light (m/s)"
                value={item.c}
                emptyBehavior="revert"
                onChange={(v) => handleChange('c', v)}
            />
            <NumberField
                label="Random Seed"
                value={item.random_seed}
                emptyBehavior="null"
                onChange={(v) => handleChange('random_seed', v)}
            />
            <NumberField
                label="ADC Bits"
                value={item.adc_bits}
                emptyBehavior="revert"
                onChange={(v) => handleChange('adc_bits', v)}
            />
            <NumberField
                label="Oversample Ratio"
                value={item.oversample_ratio}
                emptyBehavior="revert"
                onChange={(v) => handleChange('oversample_ratio', v)}
            />

            <Section title="Georeference">
                <NumberField
                    label="Origin Latitude (deg)"
                    value={item.origin.latitude}
                    emptyBehavior="revert"
                    onChange={(v) => handleChange('origin.latitude', v)}
                />
                <NumberField
                    label="Origin Longitude (deg)"
                    value={item.origin.longitude}
                    emptyBehavior="revert"
                    onChange={(v) => handleChange('origin.longitude', v)}
                />
                <NumberField
                    label="Origin Altitude (m)"
                    value={item.origin.altitude}
                    emptyBehavior="revert"
                    onChange={(v) => handleChange('origin.altitude', v)}
                />
                <FormControl fullWidth size="small">
                    <InputLabel>Coordinate System</InputLabel>
                    <Select
                        label="Coordinate System"
                        value={item.coordinateSystem.frame}
                        onChange={(e) =>
                            handleCoordinateSystemFrameChange(
                                e.target
                                    .value as GlobalParameters['coordinateSystem']['frame']
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
                            emptyBehavior="revert"
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
            <Dialog
                open={pendingRotationUnit !== null}
                onClose={() => setPendingRotationUnit(null)}
            >
                <DialogTitle>Change Rotation Angle Unit</DialogTitle>
                <DialogContent>
                    <DialogContentText>
                        Existing rotation values are present. Convert them to
                        keep the same physical orientation and rates, or keep
                        the numeric values as they are.
                    </DialogContentText>
                </DialogContent>
                <DialogActions>
                    <Button onClick={() => setPendingRotationUnit(null)}>
                        Cancel
                    </Button>
                    <Button
                        onClick={() => {
                            if (pendingRotationUnit) {
                                setRotationAngleUnit(
                                    pendingRotationUnit,
                                    false
                                );
                            }
                            setPendingRotationUnit(null);
                        }}
                    >
                        Keep Numeric Values
                    </Button>
                    <Button
                        variant="contained"
                        onClick={() => {
                            if (pendingRotationUnit) {
                                setRotationAngleUnit(pendingRotationUnit, true);
                            }
                            setPendingRotationUnit(null);
                        }}
                    >
                        Convert Existing Values
                    </Button>
                </DialogActions>
            </Dialog>
        </Box>
    );
}
