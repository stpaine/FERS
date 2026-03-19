// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import React, { useState } from 'react';
import { Box, Snackbar, Alert } from '@mui/material';
import AppRail from '@/components/AppRail';
import { ScenarioView } from '@/views/ScenarioView';
import { AssetLibraryView } from '@/views/AssetLibraryView';
import { SimulationView } from '@/views/SimulationView';
import { ResultsView } from '@/views/ResultsView';
import SettingsDialog from '@/components/SettingsDialog';
import { useScenarioStore } from '@/stores/scenarioStore';

export function MainLayout() {
    const [activeView, setActiveView] = useState('scenario');
    const [settingsOpen, setSettingsOpen] = useState(false);
    const { open: errorOpen, message: errorMessage } = useScenarioStore(
        (state) => state.errorSnackbar
    );
    const hideError = useScenarioStore((state) => state.hideError);

    return (
        <Box
            sx={{
                display: 'flex',
                height: '100vh',
                width: '100vw',
                overflow: 'hidden',
                position: 'fixed', // Ensure it stays in viewport
                top: 0,
                left: 0,
            }}
        >
            <AppRail
                activeView={activeView}
                onViewChange={setActiveView}
                onSettingsClick={() => setSettingsOpen(true)}
            />
            <Box
                component="main"
                sx={{
                    flexGrow: 1,
                    minWidth: 0, // Allow shrinking below content size
                    height: '100%',
                    overflow: 'hidden', // Prevent overflow
                    position: 'relative',
                }}
            >
                {/* Render all views but only display the active one */}
                <Box
                    sx={{
                        display: activeView === 'scenario' ? 'flex' : 'none',
                        height: '100%',
                        width: '100%',
                    }}
                >
                    <ScenarioView />
                </Box>
                <Box
                    sx={{
                        display: activeView === 'assets' ? 'block' : 'none',
                        height: '100%',
                        width: '100%',
                    }}
                >
                    <AssetLibraryView />
                </Box>
                <Box
                    sx={{
                        display: activeView === 'simulation' ? 'block' : 'none',
                        height: '100%',
                        width: '100%',
                    }}
                >
                    <SimulationView />
                </Box>
                <Box
                    sx={{
                        display: activeView === 'results' ? 'block' : 'none',
                        height: '100%',
                        width: '100%',
                    }}
                >
                    <ResultsView />
                </Box>
            </Box>
            <SettingsDialog
                open={settingsOpen}
                onClose={() => setSettingsOpen(false)}
            />
            <Snackbar
                open={errorOpen}
                autoHideDuration={6000}
                onClose={hideError}
                anchorOrigin={{ vertical: 'bottom', horizontal: 'center' }}
            >
                <Alert
                    onClose={hideError}
                    severity="error"
                    variant="filled"
                    sx={{ width: '100%' }}
                >
                    {errorMessage}
                </Alert>
            </Snackbar>
        </Box>
    );
}
