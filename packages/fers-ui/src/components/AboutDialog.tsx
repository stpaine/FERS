// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import GitHubIcon from '@mui/icons-material/GitHub';
import {
    Box,
    Button,
    Chip,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Divider,
    Typography,
} from '@mui/material';
import { getVersion } from '@tauri-apps/api/app';
import { openUrl } from '@tauri-apps/plugin-opener';
import { useEffect, useState } from 'react';

interface AboutDialogProps {
    open: boolean;
    onClose: () => void;
    onLicensesClick: () => void;
}

const GITHUB_URL = 'https://github.com/stpaine/FERS';

export default function AboutDialog({
    open,
    onClose,
    onLicensesClick,
}: AboutDialogProps) {
    const [version, setVersion] = useState('...');

    useEffect(() => {
        if (!open) return;
        getVersion()
            .then(setVersion)
            .catch(() => setVersion('unknown'));
    }, [open]);

    const handleGitHub = () => {
        openUrl(GITHUB_URL).catch(console.error);
    };

    return (
        <Dialog open={open} onClose={onClose} maxWidth="xs" fullWidth>
            <DialogTitle>About FERS</DialogTitle>
            <DialogContent>
                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'column',
                        alignItems: 'center',
                        gap: 1,
                        pt: 1,
                        pb: 2,
                        textAlign: 'center',
                    }}
                >
                    <Typography variant="h5" fontWeight="bold">
                        FERS
                    </Typography>
                    <Typography variant="subtitle2" color="text.secondary">
                        Flexible Extensible Radar Simulator
                    </Typography>
                    <Typography variant="body2">Version {version}</Typography>
                    <Chip
                        label="GPL-2.0-only"
                        size="small"
                        variant="outlined"
                        sx={{ mt: 0.5 }}
                    />
                    <Divider flexItem sx={{ my: 1 }} />
                    <Typography variant="body2" color="text.secondary">
                        Copyright © 2006–2008 Marc Brooker and Michael Inggs.
                    </Typography>
                    <Typography variant="body2" color="text.secondary">
                        Copyright © 2008–present FERS Contributors.
                    </Typography>
                    <Button
                        startIcon={<GitHubIcon />}
                        onClick={handleGitHub}
                        size="small"
                        sx={{ mt: 1 }}
                    >
                        View on GitHub
                    </Button>
                </Box>
            </DialogContent>
            <DialogActions sx={{ justifyContent: 'space-between' }}>
                <Button onClick={onLicensesClick} size="small">
                    Third-Party Licenses
                </Button>
                <Button onClick={onClose}>Close</Button>
            </DialogActions>
        </Dialog>
    );
}
