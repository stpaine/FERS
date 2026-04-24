// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import DeleteIcon from '@mui/icons-material/Delete';
import { Box, Button, IconButton } from '@mui/material';
import { Timing, useScenarioStore } from '@/stores/scenarioStore';
import { BufferedTextField, NumberField, Section } from './InspectorControls';

interface TimingInspectorProps {
    item: Timing;
}

export function TimingInspector({ item }: TimingInspectorProps) {
    const { updateItem, addNoiseEntry, removeNoiseEntry } =
        useScenarioStore.getState();
    const handleChange = (path: string, value: unknown) =>
        updateItem(item.id, path, value);

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
            <NumberField
                label="Frequency (Hz)"
                value={item.frequency}
                emptyBehavior="revert"
                onChange={(v) => handleChange('frequency', v)}
            />
            <NumberField
                label="Frequency Offset (Hz)"
                value={item.freqOffset}
                emptyBehavior="null"
                onChange={(v) => handleChange('freqOffset', v)}
            />
            <NumberField
                label="Random Freq. Offset Stdev (Hz)"
                value={item.randomFreqOffsetStdev}
                emptyBehavior="null"
                onChange={(v) => handleChange('randomFreqOffsetStdev', v)}
            />
            <NumberField
                label="Phase Offset (rad)"
                value={item.phaseOffset}
                emptyBehavior="null"
                onChange={(v) => handleChange('phaseOffset', v)}
            />
            <NumberField
                label="Random Phase Offset Stdev (rad)"
                value={item.randomPhaseOffsetStdev}
                emptyBehavior="null"
                onChange={(v) => handleChange('randomPhaseOffsetStdev', v)}
            />
            <Section title="Noise Entries">
                {item.noiseEntries.map((entry, index) => (
                    <Box
                        key={entry.id}
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
                        <Box sx={{ flexGrow: 1, display: 'flex', gap: 1 }}>
                            <NumberField
                                label="Alpha"
                                value={entry.alpha}
                                emptyBehavior="revert"
                                onChange={(v) =>
                                    handleChange(
                                        `noiseEntries.${index}.alpha`,
                                        v
                                    )
                                }
                            />
                            <NumberField
                                label="Weight"
                                value={entry.weight}
                                emptyBehavior="revert"
                                onChange={(v) =>
                                    handleChange(
                                        `noiseEntries.${index}.weight`,
                                        v
                                    )
                                }
                            />
                        </Box>
                        <IconButton
                            size="small"
                            onClick={() => removeNoiseEntry(item.id, entry.id)}
                        >
                            <DeleteIcon fontSize="small" />
                        </IconButton>
                    </Box>
                ))}
                <Button onClick={() => addNoiseEntry(item.id)} size="small">
                    Add Noise Entry
                </Button>
            </Section>
        </Box>
    );
}
