// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { createTheme, responsiveFontSizes } from '@mui/material/styles';

/**
 * FERS "Deep Radar" Color Palette.
 * A high-contrast, engineering-focused dark theme tailored for data visualization.
 */
export const fersColors = {
    background: {
        default: '#0b0e14', // Deep Gunmetal - Main app background
        paper: '#151b24', // Lighter Gunmetal - Panels and Cards
        canvas: '#020408', // Near Pitch Black - The 3D Void
        // Glassmorphism overlay for HUD elements (labels)
        overlay: 'rgba(21, 27, 36, 0.90)',
    },
    primary: {
        main: '#38bdf8', // Sky Blue - Active elements, Buttons
        dark: '#0284c7',
        light: '#bae6fd',
        contrastText: '#0f172a',
    },
    secondary: {
        main: '#fbbf24', // Amber - Selection, Highlights
        dark: '#d97706',
        light: '#fde68a',
        contrastText: '#1e1b4b',
    },
    text: {
        primary: '#f0f6fc', // Off-white for high readability
        secondary: '#8b949e', // Muted blue-grey for labels
        disabled: '#484f58',
    },
    // 3D Scene Specific Semantics
    platform: {
        default: '#38bdf8', // Cyan (Primary) - Standard Platform color
        selected: '#fbbf24', // Amber (Secondary) - Selected Platform glow
        emission: '#fbbf24', // Self-illumination color when selected
    },
    physics: {
        velocity: '#4ade80', // Neon Green - Vector forward
        boresight: '#facc15', // Bright Yellow - Where the antenna points
        motionPath: '#c084fc', // Lavender/Purple - Distinct from RF links
        rcs: '#f97316', // Orange - RCS Sphere wireframe
    },
    link: {
        monostatic: {
            strong: '#22c55e', // Emerald - Strong Signal / Lock
            weak: '#ef4444', // Red - Sub-noise / Ghost signal
        },
        illuminator: '#eab308', // Yellow - Transmit Energy
        scattered: '#06b6d4', // Cyan - Scattered Return
        direct: '#d946ef', // Fuchsia - Interference / Direct Path
    },
};

// Create the Material UI theme instance
let theme = createTheme({
    palette: {
        mode: 'dark',
        primary: {
            main: fersColors.primary.main,
            dark: fersColors.primary.dark,
            light: fersColors.primary.light,
            contrastText: fersColors.primary.contrastText,
        },
        secondary: {
            main: fersColors.secondary.main,
            dark: fersColors.secondary.dark,
            light: fersColors.secondary.light,
            contrastText: fersColors.secondary.contrastText,
        },
        background: {
            default: fersColors.background.default,
            paper: fersColors.background.paper,
        },
        text: {
            primary: fersColors.text.primary,
            secondary: fersColors.text.secondary,
            disabled: fersColors.text.disabled,
        },
        divider: '#30363d', // Subtle border color
        action: {
            hover: 'rgba(56, 189, 248, 0.08)', // Tinted with primary
            selected: 'rgba(56, 189, 248, 0.16)',
        },
    },
    typography: {
        fontFamily: ['Inter', 'Roboto', 'sans-serif'].join(','),
        h6: {
            fontWeight: 600,
            letterSpacing: '0.02em',
        },
        overline: {
            fontWeight: 700,
            letterSpacing: '0.1em',
        },
    },
    shape: {
        borderRadius: 8, // More modern, slightly rounded corners
    },
    components: {
        MuiCssBaseline: {
            styleOverrides: {
                body: {
                    scrollbarColor: '#30363d #0b0e14',
                    '&::-webkit-scrollbar, & *::-webkit-scrollbar': {
                        backgroundColor: '#0b0e14',
                        width: '8px',
                        height: '8px',
                    },
                    '&::-webkit-scrollbar-thumb, & *::-webkit-scrollbar-thumb':
                        {
                            borderRadius: 8,
                            backgroundColor: '#30363d',
                            minHeight: 24,
                            border: '2px solid #0b0e14',
                        },
                    '&::-webkit-scrollbar-thumb:focus, & *::-webkit-scrollbar-thumb:focus':
                        {
                            backgroundColor: '#8b949e',
                        },
                },
            },
        },
        MuiPaper: {
            styleOverrides: {
                root: {
                    backgroundImage: 'none', // Remove the default MUI overlay lighten effect
                    border: '1px solid #30363d', // Add thin borders for definition
                },
            },
        },
        MuiButton: {
            styleOverrides: {
                root: {
                    textTransform: 'none', // More modern text style
                    fontWeight: 600,
                },
            },
        },
    },
});

theme = responsiveFontSizes(theme);

export default theme;
