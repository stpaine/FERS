// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import LibraryAddIcon from '@mui/icons-material/LibraryAdd';
import {
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogContentText,
    DialogTitle,
    Divider,
    TextField,
    Typography,
} from '@mui/material';
import { useState } from 'react';
import { useAssetLibraryStore } from '@/stores/assetLibrary';
import {
    findComponentInStore,
    findItemInStore,
    useScenarioStore,
} from '@/stores/scenarioStore';
import { assertNever } from '@/utils/typeUtils';
import { AntennaInspector } from './inspectors/AntennaInspector';
import { GlobalParametersInspector } from './inspectors/GlobalParametersInspector';
import { PlatformComponentInspector } from './inspectors/PlatformComponentInspector';
import { PlatformInspector } from './inspectors/PlatformInspector';
import { TimingInspector } from './inspectors/TimingInspector';
import { WaveformInspector } from './inspectors/WaveformInspector';

function normalizeAssetName(name: string): string {
    return name.trim().toLocaleLowerCase();
}

function SaveToAssetLibraryButton({
    itemId,
    itemName,
}: {
    itemId: string;
    itemName: string;
}) {
    const [isSaving, setIsSaving] = useState(false);
    const [isDuplicateDialogOpen, setDuplicateDialogOpen] = useState(false);
    const [proposedName, setProposedName] = useState(itemName);
    const templates = useAssetLibraryStore((state) => state.templates);
    const hasLoaded = useAssetLibraryStore((state) => state.hasLoaded);
    const loadCatalog = useAssetLibraryStore((state) => state.loadCatalog);
    const saveScenarioItem = useAssetLibraryStore(
        (state) => state.saveScenarioItem
    );
    const showSuccess = useScenarioStore((state) => state.showSuccess);
    const showError = useScenarioStore((state) => state.showError);

    const handleSaveTemplate = async (nameOverride?: string) => {
        setIsSaving(true);
        try {
            const template = await saveScenarioItem(itemId, nameOverride);
            if (!template) {
                showError(
                    'Selected item cannot be saved as an asset template.'
                );
                return;
            }
            showSuccess(`${template.name} saved to Asset Library.`);
        } catch (error) {
            const message =
                error instanceof Error ? error.message : String(error);
            showError(`Save to Asset Library failed: ${message}`);
        } finally {
            setIsSaving(false);
        }
    };

    const handleSaveClick = async () => {
        if (!hasLoaded) {
            await loadCatalog();
        }

        const currentTemplates = useAssetLibraryStore.getState().templates;
        const duplicateExists = currentTemplates.some(
            (template) =>
                normalizeAssetName(template.name) ===
                normalizeAssetName(itemName)
        );

        if (duplicateExists) {
            setProposedName(itemName);
            setDuplicateDialogOpen(true);
            return;
        }

        await handleSaveTemplate();
    };

    const handleKeepDuplicateName = async () => {
        setDuplicateDialogOpen(false);
        await handleSaveTemplate(itemName);
    };

    const handleSaveWithNewName = async () => {
        const trimmedName = proposedName.trim();
        if (!trimmedName) {
            return;
        }

        setDuplicateDialogOpen(false);
        await handleSaveTemplate(trimmedName);
    };

    const normalizedItemName = normalizeAssetName(itemName);
    const normalizedProposedName = normalizeAssetName(proposedName);
    const proposedNameConflicts =
        normalizedProposedName.length > 0 &&
        templates.some(
            (template) =>
                normalizeAssetName(template.name) === normalizedProposedName
        );
    const canSaveWithNewName =
        normalizedProposedName.length > 0 &&
        normalizedProposedName !== normalizedItemName &&
        !proposedNameConflicts;

    return (
        <>
            <Button
                size="small"
                variant="outlined"
                startIcon={<LibraryAddIcon />}
                onClick={() => void handleSaveClick()}
                disabled={isSaving}
                sx={{ alignSelf: 'flex-start' }}
            >
                Save to Asset Library
            </Button>
            <Dialog
                open={isDuplicateDialogOpen}
                onClose={() => setDuplicateDialogOpen(false)}
            >
                <DialogTitle>Asset Name Already Exists</DialogTitle>
                <DialogContent>
                    <DialogContentText sx={{ mb: 2 }}>
                        Another saved asset already uses this name. Enter a new
                        name, or keep the current name.
                    </DialogContentText>
                    <TextField
                        autoFocus
                        fullWidth
                        label="Asset name"
                        value={proposedName}
                        onChange={(event) =>
                            setProposedName(event.target.value)
                        }
                        error={
                            proposedNameConflicts &&
                            normalizedProposedName !== normalizedItemName
                        }
                        helperText={
                            proposedNameConflicts &&
                            normalizedProposedName !== normalizedItemName
                                ? 'This name is already in use.'
                                : 'Choose a unique name for the saved asset.'
                        }
                    />
                </DialogContent>
                <DialogActions>
                    <Button onClick={() => setDuplicateDialogOpen(false)}>
                        Cancel
                    </Button>
                    <Button onClick={() => void handleKeepDuplicateName()}>
                        Keep Name
                    </Button>
                    <Button
                        variant="contained"
                        onClick={() => void handleSaveWithNewName()}
                        disabled={!canSaveWithNewName}
                    >
                        Save New Name
                    </Button>
                </DialogActions>
            </Dialog>
        </>
    );
}

function InspectorContent() {
    const selectedItemId = useScenarioStore((state) => state.selectedItemId);
    const selectedComponentId = useScenarioStore(
        (state) => state.selectedComponentId
    );
    const selectedItem = useScenarioStore((state) =>
        findItemInStore(state, selectedItemId)
    );
    const selectedComponent = useScenarioStore(
        (state) =>
            findComponentInStore(state, selectedComponentId)?.component ?? null
    );

    if (!selectedItem) {
        return (
            <Typography color="text.secondary">
                Select an item to see its properties.
            </Typography>
        );
    }

    const renderInspector = () => {
        if (selectedComponent && selectedItem.type === 'Platform') {
            const componentIndex = selectedItem.components.findIndex(
                (c) => c.id === selectedComponent.id
            );
            if (componentIndex === -1) {
                return (
                    <Typography color="text.secondary">
                        Select an item to see its properties.
                    </Typography>
                );
            }
            return (
                <PlatformComponentInspector
                    component={selectedComponent}
                    platformId={selectedItem.id}
                    index={componentIndex}
                />
            );
        }
        switch (selectedItem.type) {
            case 'GlobalParameters':
                return <GlobalParametersInspector item={selectedItem} />;
            case 'Waveform':
                return <WaveformInspector item={selectedItem} />;
            case 'Timing':
                return <TimingInspector item={selectedItem} />;
            case 'Antenna':
                return <AntennaInspector item={selectedItem} />;
            case 'Platform':
                return (
                    <PlatformInspector
                        item={selectedItem}
                        selectedComponentId={selectedComponentId}
                    />
                );
            default:
                return assertNever(selectedItem);
        }
    };

    const canSaveToAssetLibrary =
        selectedItem.type !== 'GlobalParameters' && !selectedComponent;

    return (
        <Box>
            <Typography variant="overline" color="text.secondary">
                {selectedItem.type}
            </Typography>
            <Divider sx={{ my: 1 }} />
            {canSaveToAssetLibrary && 'name' in selectedItem && (
                <Box sx={{ mb: 2 }}>
                    <SaveToAssetLibraryButton
                        itemId={selectedItem.id}
                        itemName={selectedItem.name}
                    />
                </Box>
            )}
            {renderInspector()}
        </Box>
    );
}

export default function PropertyInspector() {
    return (
        <Box
            sx={{
                height: '100%',
                display: 'flex',
                flexDirection: 'column',
            }}
        >
            <Box sx={{ flexShrink: 0, px: 2, pt: 2, pb: 1 }}>
                <Typography variant="h6">Properties</Typography>
            </Box>
            <Divider sx={{ mx: 2 }} />

            <Box
                sx={{
                    flexGrow: 1,
                    overflowY: 'auto',
                    minHeight: 0,
                    p: 1.5,
                }}
            >
                <InspectorContent />
            </Box>
        </Box>
    );
}
