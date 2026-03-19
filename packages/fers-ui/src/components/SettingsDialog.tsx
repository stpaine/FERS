// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import {
    Dialog,
    DialogTitle,
    DialogContent,
    DialogActions,
    Button,
    Typography,
    Box,
} from '@mui/material';
import { useScenarioStore } from '@/stores/scenarioStore';
import { NumberField } from './inspectors/InspectorControls';

interface SettingsDialogProps {
    open: boolean;
    onClose: () => void;
}

export default function SettingsDialog({ open, onClose }: SettingsDialogProps) {
    const { targetPlaybackDuration, setTargetPlaybackDuration } =
        useScenarioStore();
    return (
        <Dialog open={open} onClose={onClose} maxWidth="xs" fullWidth>
            <DialogTitle>Application Settings</DialogTitle>
            <DialogContent>
                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'column',
                        gap: 2,
                        pt: 1,
                    }}
                >
                    <Typography>
                        Global application settings. Scenario parameters are
                        edited in the Property Inspector.
                    </Typography>
                    <NumberField
                        label="Target Preview Playback Duration (s)"
                        value={targetPlaybackDuration}
                        onChange={(val) => setTargetPlaybackDuration(val)}
                    />
                    <Typography variant="caption">
                        Set a fixed real-world duration for the simulation
                        preview. Leave blank for default behavior (real-time
                        playback, with a minimum of 5 seconds for short
                        simulations).
                    </Typography>
                </Box>
            </DialogContent>
            <DialogActions>
                <Button onClick={onClose}>Close</Button>
            </DialogActions>
        </Dialog>
    );
}
