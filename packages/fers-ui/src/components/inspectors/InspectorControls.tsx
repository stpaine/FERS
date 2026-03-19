// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import React from 'react';
import {
    Accordion,
    AccordionSummary,
    AccordionDetails,
    Typography,
    TextField,
    Box,
} from '@mui/material';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import { open } from '@tauri-apps/plugin-dialog';

export const Section = ({
    title,
    children,
}: {
    title: string;
    children: React.ReactNode;
}) => (
    <Accordion
        defaultExpanded
        disableGutters
        elevation={0}
        sx={{
            border: 1,
            borderColor: 'divider',
            borderRadius: 1,
            '&:before': {
                display: 'none',
            },
            overflow: 'hidden',
        }}
    >
        <AccordionSummary
            expandIcon={<ExpandMoreIcon />}
            sx={{
                bgcolor: 'action.hover',
                px: 2,
                minHeight: 48,
                '& .MuiAccordionSummary-content': {
                    my: 0,
                },
                '& .MuiAccordionSummary-expandIconWrapper': {
                    color: 'text.secondary',
                },
            }}
        >
            <Typography variant="subtitle2" sx={{ fontWeight: 600 }}>
                {title}
            </Typography>
        </AccordionSummary>
        <AccordionDetails
            sx={{
                p: 2,
                display: 'flex',
                flexDirection: 'column',
                gap: 2,
                borderTop: 1,
                borderColor: 'divider',
            }}
        >
            {children}
        </AccordionDetails>
    </Accordion>
);

export const NumberField = ({
    label,
    value,
    onChange,
}: {
    label: string;
    value: number | null;
    onChange: (val: number | null) => void;
}) => {
    const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const strValue = e.target.value;

        // If the field is cleared, represent it as null.
        if (strValue === '') {
            onChange(null);
            return;
        }

        const numValue = parseFloat(strValue);

        // Only update the parent state if the parsed value is a valid number.
        // This prevents NaN from being stored and allows the user to type
        // intermediate invalid states (e.g., "1.2.3", "-") without
        // corrupting the application state.
        if (!isNaN(numValue)) {
            onChange(numValue);
        }
    };

    return (
        <TextField
            label={label}
            type="number"
            variant="outlined"
            size="small"
            fullWidth
            value={value ?? ''}
            onChange={handleChange}
            // Allow floating point numbers in the number input's spinners.
            inputProps={{ step: 'any' }}
        />
    );
};

export const FileInput = ({
    label,
    value,
    onChange,
    filters,
}: {
    label: string;
    value?: string;
    onChange: (val: string) => void;
    filters: { name: string; extensions: string[] }[];
}) => {
    const handleOpenFile = async () => {
        try {
            const selected = await open({
                multiple: false,
                filters,
            });
            if (typeof selected === 'string') {
                onChange(selected);
            }
        } catch (err) {
            console.error('Error opening file dialog:', err);
        }
    };

    return (
        <Box>
            <Typography variant="caption" color="text.secondary">
                {label}
            </Typography>
            <TextField
                variant="outlined"
                size="small"
                fullWidth
                value={value?.split(/[/\\]/).pop() ?? 'No file selected.'}
                onClick={handleOpenFile}
                sx={{ cursor: 'pointer' }}
                slotProps={{
                    input: {
                        readOnly: true,
                        style: { cursor: 'pointer' },
                    },
                }}
            />
        </Box>
    );
};
