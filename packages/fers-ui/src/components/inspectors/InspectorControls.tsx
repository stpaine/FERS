// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import {
    Accordion,
    AccordionDetails,
    AccordionSummary,
    Box,
    TextField,
    Typography,
} from '@mui/material';
import { open } from '@tauri-apps/plugin-dialog';
import React from 'react';

export type EmptyFieldBehavior = 'revert' | 'null';

type NumberBlurResolution =
    | {
          action: 'commit';
          nextDraft: string;
          value: number | null;
      }
    | {
          action: 'revert';
          error: string;
          nextDraft: string;
      };

type TextBlurResolution =
    | {
          action: 'commit';
          nextDraft: string;
          value: string;
      }
    | {
          action: 'revert';
          error: string;
          nextDraft: string;
      };

export function formatNumberFieldValue(value: number | null): string {
    return value === null ? '' : String(value);
}

export function resolveNumberFieldBlur(
    draft: string,
    committedValue: number | null,
    emptyBehavior: EmptyFieldBehavior
): NumberBlurResolution {
    const trimmedDraft = draft.trim();

    if (trimmedDraft === '') {
        if (emptyBehavior === 'null') {
            return {
                action: 'commit',
                nextDraft: '',
                value: null,
            };
        }

        return {
            action: 'revert',
            error: 'Value required.',
            nextDraft: formatNumberFieldValue(committedValue),
        };
    }

    const parsedValue = Number(trimmedDraft);
    if (!Number.isFinite(parsedValue)) {
        return {
            action: 'revert',
            error: 'Invalid number.',
            nextDraft: formatNumberFieldValue(committedValue),
        };
    }

    return {
        action: 'commit',
        nextDraft: String(parsedValue),
        value: parsedValue,
    };
}

export function resolveTextFieldBlur(
    draft: string,
    committedValue: string,
    allowEmpty: boolean
): TextBlurResolution {
    if (!allowEmpty && draft.trim() === '') {
        return {
            action: 'revert',
            error: 'Value required.',
            nextDraft: committedValue,
        };
    }

    return {
        action: 'commit',
        nextDraft: draft,
        value: draft,
    };
}

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
    emptyBehavior = 'revert',
    helperText,
    externalError = false,
}: {
    label: string;
    value: number | null;
    onChange: (val: number | null) => void;
    emptyBehavior?: EmptyFieldBehavior;
    helperText?: string;
    externalError?: boolean;
}) => {
    const [draft, setDraft] = React.useState(() =>
        formatNumberFieldValue(value)
    );
    const [isFocused, setIsFocused] = React.useState(false);
    const [error, setError] = React.useState<string | null>(null);

    React.useEffect(() => {
        if (!isFocused) {
            setDraft(formatNumberFieldValue(value));
        }
    }, [isFocused, value]);

    const handleBlur = () => {
        setIsFocused(false);

        const resolution = resolveNumberFieldBlur(draft, value, emptyBehavior);
        setDraft(resolution.nextDraft);

        if (resolution.action === 'commit') {
            setError(null);
            if (resolution.value !== value) {
                onChange(resolution.value);
            }
            return;
        }

        setError(resolution.error);
    };

    return (
        <TextField
            label={label}
            type="text"
            variant="outlined"
            size="small"
            fullWidth
            value={draft}
            error={error !== null || externalError}
            helperText={error ?? helperText}
            onFocus={() => {
                setIsFocused(true);
                setError(null);
            }}
            onBlur={handleBlur}
            onChange={(e) => {
                setDraft(e.target.value);
                if (error) {
                    setError(null);
                }
            }}
            onKeyDown={(e) => {
                if (e.key === 'Enter') {
                    e.currentTarget.blur();
                } else if (e.key === 'Escape') {
                    setDraft(formatNumberFieldValue(value));
                    setError(null);
                    e.currentTarget.blur();
                }
            }}
            inputProps={{
                inputMode: 'decimal',
            }}
        />
    );
};

export const BufferedTextField = ({
    label,
    value,
    onChange,
    allowEmpty = true,
    ...textFieldProps
}: Omit<
    React.ComponentProps<typeof TextField>,
    'label' | 'value' | 'onChange'
> & {
    label: string;
    value: string;
    onChange: (val: string) => void;
    allowEmpty?: boolean;
}) => {
    const [draft, setDraft] = React.useState(value);
    const [isFocused, setIsFocused] = React.useState(false);
    const [error, setError] = React.useState<string | null>(null);

    React.useEffect(() => {
        if (!isFocused) {
            setDraft(value);
        }
    }, [isFocused, value]);

    const handleBlur = () => {
        setIsFocused(false);

        const resolution = resolveTextFieldBlur(draft, value, allowEmpty);
        setDraft(resolution.nextDraft);

        if (resolution.action === 'commit') {
            setError(null);
            if (resolution.value !== value) {
                onChange(resolution.value);
            }
            return;
        }

        setError(resolution.error);
    };

    return (
        <TextField
            {...textFieldProps}
            label={label}
            value={draft}
            error={error !== null || textFieldProps.error === true}
            helperText={error ?? textFieldProps.helperText}
            onFocus={(e) => {
                setIsFocused(true);
                setError(null);
                textFieldProps.onFocus?.(e);
            }}
            onBlur={(e) => {
                handleBlur();
                textFieldProps.onBlur?.(e);
            }}
            onChange={(e) => {
                setDraft(e.target.value);
                if (error) {
                    setError(null);
                }
            }}
            onKeyDown={(e) => {
                if (e.key === 'Enter') {
                    e.currentTarget.blur();
                } else if (e.key === 'Escape') {
                    setDraft(value);
                    setError(null);
                    e.currentTarget.blur();
                }
                textFieldProps.onKeyDown?.(e);
            }}
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
