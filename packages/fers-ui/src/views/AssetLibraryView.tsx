// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import DeleteIcon from '@mui/icons-material/Delete';
import EditIcon from '@mui/icons-material/Edit';
import FileDownloadIcon from '@mui/icons-material/FileDownload';
import FileUploadIcon from '@mui/icons-material/FileUpload';
import LibraryAddIcon from '@mui/icons-material/LibraryAdd';
import SaveIcon from '@mui/icons-material/Save';
import {
    Alert,
    Box,
    Button,
    Checkbox,
    Chip,
    Divider,
    IconButton,
    Paper,
    Stack,
    TextField,
    Tooltip,
    Typography,
} from '@mui/material';
import { open, save } from '@tauri-apps/plugin-dialog';
import { readTextFile, writeTextFile } from '@tauri-apps/plugin-fs';
import { useEffect, useMemo, useState } from 'react';
import ConfirmDialog from '@/components/ConfirmDialog';
import {
    createAssetLibraryFile,
    parseAssetTemplates,
    useAssetLibraryStore,
} from '@/stores/assetLibrary';
import type { AssetLibraryTemplate } from '@/stores/assetTemplates';
import { useScenarioStore } from '@/stores/scenarioStore';

const assetKindLabels: Record<AssetLibraryTemplate['kind'], string> = {
    waveform: 'Waveform',
    timing: 'Timing',
    antenna: 'Antenna',
    platform: 'Platform',
};

const assetKindOrder: AssetLibraryTemplate['kind'][] = [
    'waveform',
    'timing',
    'antenna',
    'platform',
];

function formatTimestamp(timestamp: string): string {
    const date = new Date(timestamp);
    if (Number.isNaN(date.getTime())) {
        return timestamp;
    }

    return date.toLocaleString();
}

function describeTemplate(template: AssetLibraryTemplate): string {
    switch (template.kind) {
        case 'waveform':
            return `${template.payload.waveformType}, ${template.payload.carrier_frequency} Hz`;
        case 'timing':
            return `${template.payload.frequency} Hz timing`;
        case 'antenna':
            return `${template.payload.pattern} antenna`;
        case 'platform': {
            const dependencyCount =
                template.dependencies.waveforms.length +
                template.dependencies.timings.length +
                template.dependencies.antennas.length;
            return `${template.payload.components.length} components, ${dependencyCount} dependencies`;
        }
    }
}

export function AssetLibraryView() {
    const templates = useAssetLibraryStore((state) => state.templates);
    const isLoading = useAssetLibraryStore((state) => state.isLoading);
    const error = useAssetLibraryStore((state) => state.error);
    const loadCatalog = useAssetLibraryStore((state) => state.loadCatalog);
    const addTemplates = useAssetLibraryStore((state) => state.addTemplates);
    const updateTemplateName = useAssetLibraryStore(
        (state) => state.updateTemplateName
    );
    const deleteTemplate = useAssetLibraryStore(
        (state) => state.deleteTemplate
    );
    const insertAssetTemplate = useScenarioStore(
        (state) => state.insertAssetTemplate
    );
    const showSuccess = useScenarioStore((state) => state.showSuccess);
    const showError = useScenarioStore((state) => state.showError);

    const [selectedIds, setSelectedIds] = useState<Set<string>>(new Set());
    const [renamingId, setRenamingId] = useState<string | null>(null);
    const [renameValue, setRenameValue] = useState('');
    const [searchQuery, setSearchQuery] = useState('');
    const [deleteCandidate, setDeleteCandidate] =
        useState<AssetLibraryTemplate | null>(null);

    useEffect(() => {
        void loadCatalog();
    }, [loadCatalog]);

    const selectedTemplates = useMemo(
        () => templates.filter((template) => selectedIds.has(template.id)),
        [selectedIds, templates]
    );

    const filteredTemplates = useMemo(() => {
        const normalizedQuery = searchQuery.trim().toLocaleLowerCase();
        if (!normalizedQuery) {
            return templates;
        }

        return templates.filter((template) =>
            template.name.toLocaleLowerCase().includes(normalizedQuery)
        );
    }, [searchQuery, templates]);

    const filteredTemplateIds = useMemo(
        () => filteredTemplates.map((template) => template.id),
        [filteredTemplates]
    );

    const allFilteredSelected =
        filteredTemplateIds.length > 0 &&
        filteredTemplateIds.every((templateId) => selectedIds.has(templateId));
    const someFilteredSelected = filteredTemplateIds.some((templateId) =>
        selectedIds.has(templateId)
    );

    const templatesByKind = useMemo(() => {
        return assetKindOrder.map((kind) => ({
            kind,
            templates: filteredTemplates.filter(
                (template) => template.kind === kind
            ),
        }));
    }, [filteredTemplates]);

    const handleToggleSelection = (templateId: string) => {
        setSelectedIds((current) => {
            const next = new Set(current);
            if (next.has(templateId)) {
                next.delete(templateId);
            } else {
                next.add(templateId);
            }
            return next;
        });
    };

    const handleSelectAllFiltered = () => {
        setSelectedIds((current) => {
            const next = new Set(current);
            filteredTemplateIds.forEach((templateId) => next.add(templateId));
            return next;
        });
    };

    const handleDeselectFiltered = () => {
        setSelectedIds((current) => {
            const next = new Set(current);
            filteredTemplateIds.forEach((templateId) =>
                next.delete(templateId)
            );
            return next;
        });
    };

    const handleImport = async () => {
        try {
            const selectedPath = await open({
                title: 'Import Asset Templates',
                multiple: false,
                filters: [
                    {
                        name: 'FERS Asset Templates',
                        extensions: ['json'],
                    },
                ],
            });

            if (typeof selectedPath !== 'string') {
                return;
            }

            const raw = await readTextFile(selectedPath);
            const importedTemplates = parseAssetTemplates(JSON.parse(raw));
            const count = await addTemplates(importedTemplates);
            showSuccess(
                `Imported ${count} asset template${count === 1 ? '' : 's'}.`
            );
        } catch (importError) {
            const message =
                importError instanceof Error
                    ? importError.message
                    : String(importError);
            showError(`Asset import failed: ${message}`);
        }
    };

    const handleExportSelected = async () => {
        if (selectedTemplates.length === 0) {
            return;
        }

        try {
            const filePath = await save({
                title: 'Export Asset Templates',
                defaultPath: 'fers-assets.fersasset.json',
                filters: [
                    {
                        name: 'FERS Asset Templates',
                        extensions: ['json'],
                    },
                ],
            });

            if (!filePath) {
                return;
            }

            await writeTextFile(
                filePath,
                JSON.stringify(
                    createAssetLibraryFile(selectedTemplates),
                    null,
                    2
                )
            );
            showSuccess(
                `Exported ${selectedTemplates.length} asset template${selectedTemplates.length === 1 ? '' : 's'}.`
            );
        } catch (exportError) {
            const message =
                exportError instanceof Error
                    ? exportError.message
                    : String(exportError);
            showError(`Asset export failed: ${message}`);
        }
    };

    const handleLoadTemplate = (template: AssetLibraryTemplate) => {
        const result = insertAssetTemplate(template);
        showSuccess(
            `${result.insertedName ?? template.name} loaded into scenario.`
        );
    };

    const handleStartRename = (template: AssetLibraryTemplate) => {
        setRenamingId(template.id);
        setRenameValue(template.name);
    };

    const handleSaveRename = async () => {
        if (!renamingId) {
            return;
        }

        try {
            await updateTemplateName(renamingId, renameValue);
            setRenamingId(null);
            setRenameValue('');
            showSuccess('Asset template renamed.');
        } catch (renameError) {
            const message =
                renameError instanceof Error
                    ? renameError.message
                    : String(renameError);
            showError(`Rename failed: ${message}`);
        }
    };

    const handleConfirmDelete = async () => {
        if (!deleteCandidate) {
            return;
        }

        try {
            await deleteTemplate(deleteCandidate.id);
            setSelectedIds((current) => {
                const next = new Set(current);
                next.delete(deleteCandidate.id);
                return next;
            });
            showSuccess('Asset template deleted.');
        } catch (deleteError) {
            const message =
                deleteError instanceof Error
                    ? deleteError.message
                    : String(deleteError);
            showError(`Delete failed: ${message}`);
        } finally {
            setDeleteCandidate(null);
        }
    };

    return (
        <Box
            sx={{
                height: '100%',
                width: '100%',
                overflow: 'auto',
                p: 3,
            }}
        >
            <Stack
                direction={{ xs: 'column', sm: 'row' }}
                spacing={2}
                sx={{ mb: 3 }}
                alignItems={{ xs: 'stretch', sm: 'center' }}
                justifyContent="space-between"
            >
                <Box>
                    <Typography variant="h4">Asset Library</Typography>
                    <Typography color="text.secondary">
                        Save reusable waveforms, timings, antennas, and
                        platforms for new scenarios.
                    </Typography>
                </Box>
                <Stack direction="row" spacing={1}>
                    <Button
                        variant="outlined"
                        startIcon={<FileUploadIcon />}
                        onClick={handleImport}
                    >
                        Import
                    </Button>
                    <Button
                        variant="contained"
                        startIcon={<FileDownloadIcon />}
                        onClick={handleExportSelected}
                        disabled={selectedTemplates.length === 0}
                    >
                        Export Selected
                    </Button>
                </Stack>
            </Stack>

            {error && (
                <Alert severity="warning" sx={{ mb: 2 }}>
                    {error}
                </Alert>
            )}

            <Stack
                direction={{ xs: 'column', sm: 'row' }}
                spacing={1}
                alignItems={{ xs: 'stretch', sm: 'center' }}
                sx={{ mb: 2 }}
            >
                <TextField
                    label="Search by name"
                    value={searchQuery}
                    onChange={(event) => setSearchQuery(event.target.value)}
                    size="small"
                    fullWidth
                />
                <Button
                    variant="outlined"
                    onClick={handleSelectAllFiltered}
                    disabled={
                        filteredTemplateIds.length === 0 || allFilteredSelected
                    }
                    sx={{ whiteSpace: 'nowrap' }}
                >
                    Select All
                </Button>
                <Button
                    variant="outlined"
                    onClick={handleDeselectFiltered}
                    disabled={!someFilteredSelected}
                    sx={{ whiteSpace: 'nowrap' }}
                >
                    Deselect
                </Button>
            </Stack>

            {isLoading ? (
                <Alert severity="info">Loading asset library.</Alert>
            ) : templates.length === 0 ? (
                <Alert severity="info">
                    Save assets from the Scenario properties panel to build your
                    library.
                </Alert>
            ) : filteredTemplates.length === 0 ? (
                <Alert severity="info">No saved assets match that name.</Alert>
            ) : (
                <Stack spacing={2}>
                    {templatesByKind.map(
                        ({ kind, templates: kindTemplates }) =>
                            kindTemplates.length > 0 ? (
                                <Box key={kind}>
                                    <Typography
                                        variant="overline"
                                        color="text.secondary"
                                    >
                                        {assetKindLabels[kind]}
                                    </Typography>
                                    <Stack spacing={1}>
                                        {kindTemplates.map((template) => (
                                            <Paper
                                                key={template.id}
                                                variant="outlined"
                                                sx={{
                                                    p: 1,
                                                    borderRadius: 1,
                                                }}
                                            >
                                                <Stack
                                                    direction={{
                                                        xs: 'column',
                                                        sm: 'row',
                                                    }}
                                                    spacing={1}
                                                    alignItems={{
                                                        xs: 'stretch',
                                                        sm: 'center',
                                                    }}
                                                >
                                                    <Checkbox
                                                        size="small"
                                                        checked={selectedIds.has(
                                                            template.id
                                                        )}
                                                        onChange={() =>
                                                            handleToggleSelection(
                                                                template.id
                                                            )
                                                        }
                                                        sx={{
                                                            p: 0.5,
                                                            alignSelf: {
                                                                xs: 'flex-start',
                                                                sm: 'center',
                                                            },
                                                        }}
                                                    />
                                                    <Box
                                                        sx={{
                                                            flex: 1,
                                                            minWidth: 0,
                                                        }}
                                                    >
                                                        <Stack
                                                            direction="row"
                                                            spacing={1}
                                                            alignItems="center"
                                                            sx={{ mb: 0.25 }}
                                                        >
                                                            <Chip
                                                                label={
                                                                    assetKindLabels[
                                                                        template
                                                                            .kind
                                                                    ]
                                                                }
                                                                size="small"
                                                                sx={{
                                                                    height: 22,
                                                                }}
                                                            />
                                                            <Typography
                                                                variant="caption"
                                                                color="text.secondary"
                                                            >
                                                                Updated{' '}
                                                                {formatTimestamp(
                                                                    template.updatedAt
                                                                )}
                                                            </Typography>
                                                        </Stack>
                                                        {renamingId ===
                                                        template.id ? (
                                                            <Stack
                                                                direction="row"
                                                                spacing={0.5}
                                                            >
                                                                <TextField
                                                                    size="small"
                                                                    value={
                                                                        renameValue
                                                                    }
                                                                    onChange={(
                                                                        event
                                                                    ) =>
                                                                        setRenameValue(
                                                                            event
                                                                                .target
                                                                                .value
                                                                        )
                                                                    }
                                                                    fullWidth
                                                                />
                                                                <IconButton
                                                                    size="small"
                                                                    aria-label="Save asset template name"
                                                                    onClick={() =>
                                                                        void handleSaveRename()
                                                                    }
                                                                >
                                                                    <SaveIcon />
                                                                </IconButton>
                                                            </Stack>
                                                        ) : (
                                                            <Typography
                                                                variant="subtitle1"
                                                                sx={{
                                                                    overflow:
                                                                        'hidden',
                                                                    textOverflow:
                                                                        'ellipsis',
                                                                }}
                                                            >
                                                                {template.name}
                                                            </Typography>
                                                        )}
                                                        <Typography
                                                            variant="body2"
                                                            color="text.secondary"
                                                            sx={{
                                                                overflow:
                                                                    'hidden',
                                                                textOverflow:
                                                                    'ellipsis',
                                                            }}
                                                        >
                                                            {describeTemplate(
                                                                template
                                                            )}
                                                        </Typography>
                                                    </Box>
                                                    <Divider
                                                        flexItem
                                                        orientation="vertical"
                                                        sx={{
                                                            display: {
                                                                xs: 'none',
                                                                sm: 'block',
                                                            },
                                                        }}
                                                    />
                                                    <Stack
                                                        direction="row"
                                                        spacing={0.5}
                                                        justifyContent="flex-end"
                                                    >
                                                        <Button
                                                            size="small"
                                                            variant="contained"
                                                            startIcon={
                                                                <LibraryAddIcon />
                                                            }
                                                            onClick={() =>
                                                                handleLoadTemplate(
                                                                    template
                                                                )
                                                            }
                                                        >
                                                            Load
                                                        </Button>
                                                        <Tooltip title="Rename">
                                                            <IconButton
                                                                size="small"
                                                                onClick={() =>
                                                                    handleStartRename(
                                                                        template
                                                                    )
                                                                }
                                                            >
                                                                <EditIcon />
                                                            </IconButton>
                                                        </Tooltip>
                                                        <Tooltip title="Delete">
                                                            <IconButton
                                                                size="small"
                                                                color="error"
                                                                onClick={() =>
                                                                    setDeleteCandidate(
                                                                        template
                                                                    )
                                                                }
                                                            >
                                                                <DeleteIcon />
                                                            </IconButton>
                                                        </Tooltip>
                                                    </Stack>
                                                </Stack>
                                            </Paper>
                                        ))}
                                    </Stack>
                                </Box>
                            ) : null
                    )}
                </Stack>
            )}

            <ConfirmDialog
                open={deleteCandidate !== null}
                title="Delete Asset Template?"
                message={`Delete ${deleteCandidate?.name ?? 'this asset template'} from the library?`}
                onConfirm={() => void handleConfirmDelete()}
                onCancel={() => setDeleteCandidate(null)}
            />
        </Box>
    );
}
