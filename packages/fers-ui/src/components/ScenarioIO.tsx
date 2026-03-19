// SPDX-License-Identifier: GPL-2-0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { IconButton, Tooltip } from '@mui/material';
import FileUploadIcon from '@mui/icons-material/FileUpload';
import FileDownloadIcon from '@mui/icons-material/FileDownload';
import { useScenarioStore } from '@/stores/scenarioStore';
import { invoke } from '@tauri-apps/api/core';
import { save, open } from '@tauri-apps/plugin-dialog';
import { writeTextFile } from '@tauri-apps/plugin-fs';
import { useState } from 'react';
import ConfirmDialog from './ConfirmDialog';

export default function ScenarioIO() {
    const loadScenario = useScenarioStore((state) => state.loadScenario);
    const isDirty = useScenarioStore((state) => state.isDirty);
    const resetScenario = useScenarioStore((state) => state.resetScenario);
    const showError = useScenarioStore((state) => state.showError);

    const [isConfirmOpen, setConfirmOpen] = useState(false);

    const handleExport = async () => {
        try {
            await useScenarioStore.getState().syncBackend();

            const xmlContent = await invoke<string>('get_scenario_as_xml');

            const filePath = await save({
                title: 'Export Scenario',
                filters: [
                    {
                        name: 'FERS XML Scenario',
                        extensions: ['xml', 'fersxml'],
                    },
                ],
            });

            if (filePath) {
                await writeTextFile(filePath, xmlContent);
                console.log('Scenario exported successfully to:', filePath);
            }
        } catch (error) {
            const errorMessage =
                error instanceof Error ? error.message : String(error);
            console.error('Failed to export scenario:', errorMessage);
            showError(`Export failed: ${errorMessage}`);
        }
    };

    const performImport = async () => {
        try {
            const selectedPath = await open({
                title: 'Import Scenario',
                multiple: false,
                filters: [
                    {
                        name: 'FERS XML Scenario',
                        extensions: ['xml', 'fersxml'],
                    },
                ],
            });

            if (typeof selectedPath === 'string') {
                // Load the XML file into the C++ core
                await invoke('load_scenario_from_xml_file', {
                    filepath: selectedPath,
                });

                // Fetch the new state as JSON from the C++ core
                const jsonState = await invoke<string>('get_scenario_as_json');
                const scenarioData = JSON.parse(jsonState);

                // Update the UI's Zustand store with the new state after resetting the current state
                resetScenario();
                loadScenario(scenarioData);

                console.log(
                    'Scenario imported and synchronized successfully from:',
                    selectedPath
                );
            }
        } catch (error) {
            const errorMessage =
                error instanceof Error ? error.message : String(error);
            console.error('Failed to import scenario:', errorMessage);
            showError(`Import failed: ${errorMessage}`);
        }
    };

    const handleImport = () => {
        if (isDirty) {
            setConfirmOpen(true);
        } else {
            void performImport();
        }
    };

    const handleConfirmImport = () => {
        setConfirmOpen(false);
        void performImport();
    };

    const handleCancelImport = () => {
        setConfirmOpen(false);
    };

    return (
        <>
            <Tooltip title="Import Scenario (XML)">
                <IconButton size="small" onClick={handleImport}>
                    <FileUploadIcon fontSize="inherit" />
                </IconButton>
            </Tooltip>
            <Tooltip title="Export Scenario (XML)">
                <IconButton size="small" onClick={handleExport}>
                    <FileDownloadIcon fontSize="inherit" />
                </IconButton>
            </Tooltip>
            <ConfirmDialog
                open={isConfirmOpen}
                onConfirm={handleConfirmImport}
                onCancel={handleCancelImport}
                title="Overwrite Current Scenario?"
                message="Importing a new scenario will discard all unsaved changes. Are you sure you want to proceed?"
            />
        </>
    );
}
