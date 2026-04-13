// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { Alert, Box, Snackbar } from '@mui/material';
import { listen } from '@tauri-apps/api/event';
import React, { useEffect, useState } from 'react';
import AppRail from '@/components/AppRail';
import RawLogDrawer from '@/components/RawLogDrawer';
import SettingsDialog from '@/components/SettingsDialog';
import { type FersLogEntry, useFersLogStore } from '@/stores/fersLogStore';
import { useScenarioStore } from '@/stores/scenarioStore';
import { AssetLibraryView } from '@/views/AssetLibraryView';
import { ResultsView } from '@/views/ResultsView';
import { ScenarioView } from '@/views/ScenarioView';
import { SimulationView } from '@/views/SimulationView';

export function MainLayout() {
    const [activeView, setActiveView] = useState('scenario');
    const [settingsOpen, setSettingsOpen] = useState(false);
    const logOpen = useFersLogStore((state) => state.isOpen);
    const appendLog = useFersLogStore((state) => state.appendLog);
    const { open, message, severity } = useScenarioStore(
        (state) => state.notificationSnackbar
    );
    const hideNotification = useScenarioStore(
        (state) => state.hideNotification
    );
    const advanceNotification = useScenarioStore(
        (state) => state.advanceNotification
    );

    useEffect(() => {
        let active = true;
        let unlistenLog: (() => void) | undefined;

        listen<FersLogEntry>('fers-log', (event) => {
            appendLog(event.payload);
        }).then((unlisten) => {
            if (active) {
                unlistenLog = unlisten;
            } else {
                unlisten();
            }
        });

        return () => {
            active = false;
            unlistenLog?.();
        };
    }, [appendLog]);

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
            {logOpen && <RawLogDrawer />}
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
                    <ScenarioView isActive={activeView === 'scenario'} />
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
                open={open}
                autoHideDuration={6000}
                onClose={hideNotification}
                anchorOrigin={{ vertical: 'bottom', horizontal: 'center' }}
                TransitionProps={{ onExited: advanceNotification }}
            >
                <Alert
                    onClose={hideNotification}
                    severity={severity}
                    variant="filled"
                    sx={{ width: '100%' }}
                >
                    {message}
                </Alert>
            </Snackbar>
        </Box>
    );
}
