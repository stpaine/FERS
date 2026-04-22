// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import {
    Accordion,
    AccordionDetails,
    AccordionSummary,
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Tab,
    Tabs,
    Typography,
} from '@mui/material';
import { useState } from 'react';
import catch2License from '../../../../THIRD_PARTY_LICENSES/catch2-LICENSE.txt?raw';
import geoLibLicense from '../../../../THIRD_PARTY_LICENSES/GeographicLib-LICENSE.txt?raw';
import highFiveLicense from '../../../../THIRD_PARTY_LICENSES/HighFive-LICENSE.txt?raw';
import jsLicenses from '../../../../THIRD_PARTY_LICENSES/js-licenses.txt?raw';
import hdf5License from '../../../../THIRD_PARTY_LICENSES/libhdf5-LICENSE.txt?raw';
import libxml2License from '../../../../THIRD_PARTY_LICENSES/libxml2-LICENSE.txt?raw';
import nlohmannLicense from '../../../../THIRD_PARTY_LICENSES/nlohmann-json-LICENSE.txt?raw';
import rustLicenses from '../../../../THIRD_PARTY_LICENSES/rust-licenses.html?raw';

interface LicensesDialogProps {
    open: boolean;
    onClose: () => void;
}

const cppLibraries = [
    { name: 'Catch2', content: catch2License },
    { name: 'GeographicLib', content: geoLibLicense },
    { name: 'HDF5', content: hdf5License },
    { name: 'HighFive', content: highFiveLicense },
    { name: 'libxml2', content: libxml2License },
    { name: 'nlohmann/json', content: nlohmannLicense },
];

const licensePreStyle: React.CSSProperties = {
    margin: 0,
    whiteSpace: 'pre-wrap',
    wordBreak: 'break-word',
    fontSize: '0.75rem',
};

export default function LicensesDialog({ open, onClose }: LicensesDialogProps) {
    const [tab, setTab] = useState(0);

    return (
        <Dialog
            open={open}
            onClose={onClose}
            maxWidth="md"
            fullWidth
            PaperProps={{ sx: { height: '80vh' } }}
        >
            <DialogTitle>Third-Party Licenses</DialogTitle>
            <DialogContent
                sx={{
                    display: 'flex',
                    flexDirection: 'column',
                    p: 0,
                    overflow: 'hidden',
                }}
            >
                <Tabs
                    value={tab}
                    onChange={(_e, v: number) => setTab(v)}
                    sx={{ borderBottom: 1, borderColor: 'divider', px: 2 }}
                >
                    <Tab label="JavaScript" />
                    <Tab label="Rust" />
                    <Tab label="C++" />
                </Tabs>
                {tab === 0 && (
                    <Box sx={{ flex: 1, overflow: 'auto', p: 1 }}>
                        <pre style={licensePreStyle}>{jsLicenses}</pre>
                    </Box>
                )}
                {tab === 1 && (
                    <Box sx={{ flex: 1, overflow: 'hidden' }}>
                        <iframe
                            srcDoc={rustLicenses}
                            style={{
                                width: '100%',
                                height: '100%',
                                border: 'none',
                            }}
                            sandbox="allow-same-origin"
                            title="Rust third-party licenses"
                        />
                    </Box>
                )}
                {tab === 2 && (
                    <Box sx={{ flex: 1, overflow: 'auto' }}>
                        {cppLibraries.map(({ name, content }) => (
                            <Accordion
                                key={name}
                                disableGutters
                                elevation={0}
                                square
                            >
                                <AccordionSummary
                                    expandIcon={<ExpandMoreIcon />}
                                >
                                    <Typography variant="body2">
                                        {name}
                                    </Typography>
                                </AccordionSummary>
                                <AccordionDetails>
                                    <pre style={licensePreStyle}>{content}</pre>
                                </AccordionDetails>
                            </Accordion>
                        ))}
                    </Box>
                )}
            </DialogContent>
            <DialogActions>
                <Button onClick={onClose}>Close</Button>
            </DialogActions>
        </Dialog>
    );
}
