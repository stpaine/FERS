// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import AddIcon from '@mui/icons-material/Add';
import AdjustIcon from '@mui/icons-material/Adjust';
import ChevronRightIcon from '@mui/icons-material/ChevronRight';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import FlightIcon from '@mui/icons-material/Flight';
import PodcastsIcon from '@mui/icons-material/Podcasts';
import PublicIcon from '@mui/icons-material/Public';
import RemoveIcon from '@mui/icons-material/Remove';
import RssFeedIcon from '@mui/icons-material/RssFeed';
import SensorsIcon from '@mui/icons-material/Sensors';
import SettingsInputAntennaIcon from '@mui/icons-material/SettingsInputAntenna';
import TimerIcon from '@mui/icons-material/Timer';
import WavesIcon from '@mui/icons-material/Waves';
import { Box, Divider, IconButton, Tooltip, Typography } from '@mui/material';
import { SimpleTreeView, TreeItem } from '@mui/x-tree-view';
import React from 'react';
import { useShallow } from 'zustand/react/shallow';
import {
    Platform,
    PlatformComponent,
    useScenarioStore,
} from '@/stores/scenarioStore';
import ScenarioIO from './ScenarioIO';

const SectionHeader = ({
    title,
    onAdd,
}: {
    title: string;
    onAdd: () => void;
}) => (
    <Box
        sx={{
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            width: '100%',
            minWidth: 0, // Allow the component to shrink within its flex container
        }}
    >
        <Typography
            variant="overline"
            sx={{
                color: 'text.secondary',
                whiteSpace: 'nowrap',
                overflow: 'hidden',
                textOverflow: 'ellipsis',
                pr: 1, // Padding between text and add button
            }}
        >
            {title}
        </Typography>
        <Tooltip title={`Add ${title.slice(0, -1)}`}>
            <IconButton
                size="small"
                onClick={(e) => {
                    e.stopPropagation();
                    onAdd();
                }}
            >
                <AddIcon fontSize="inherit" />
            </IconButton>
        </Tooltip>
    </Box>
);

const getPlatformIcon = (platform: Platform) => {
    const types = platform.components.map((c) => c.type);

    if (types.includes('monostatic')) {
        return <SensorsIcon sx={{ mr: 1 }} fontSize="small" />;
    }
    if (types.includes('transmitter')) {
        return <PodcastsIcon sx={{ mr: 1 }} fontSize="small" />;
    }
    if (types.includes('receiver')) {
        return <RssFeedIcon sx={{ mr: 1 }} fontSize="small" />;
    }
    if (types.includes('target')) {
        return <AdjustIcon sx={{ mr: 1 }} fontSize="small" />;
    }
    return <FlightIcon sx={{ mr: 1 }} fontSize="small" />;
};

const getComponentIcon = (component: PlatformComponent) => {
    switch (component.type) {
        case 'monostatic':
            return <SensorsIcon sx={{ mr: 1 }} fontSize="small" />;
        case 'transmitter':
            return <PodcastsIcon sx={{ mr: 1 }} fontSize="small" />;
        case 'receiver':
            return <RssFeedIcon sx={{ mr: 1 }} fontSize="small" />;
        case 'target':
            return <AdjustIcon sx={{ mr: 1 }} fontSize="small" />;
        default:
            return <FlightIcon sx={{ mr: 1 }} fontSize="small" />;
    }
};

export default function SceneTree() {
    const {
        globalParameters,
        waveforms,
        timings,
        antennas,
        platforms,
        selectedItemId,
        selectedComponentId,
        selectItem,
        addWaveform,
        addTiming,
        addAntenna,
        addPlatform,
        removeItem,
        removePlatformComponent,
    } = useScenarioStore(
        useShallow((state) => ({
            globalParameters: state.globalParameters,
            waveforms: state.waveforms,
            timings: state.timings,
            antennas: state.antennas,
            platforms: state.platforms,
            selectedItemId: state.selectedItemId,
            selectedComponentId: state.selectedComponentId,
            selectItem: state.selectItem,
            addWaveform: state.addWaveform,
            addTiming: state.addTiming,
            addAntenna: state.addAntenna,
            addPlatform: state.addPlatform,
            removeItem: state.removeItem,
            removePlatformComponent: state.removePlatformComponent,
        }))
    );

    const handleSelect = (
        _event: React.SyntheticEvent | null,
        nodeId: string | null
    ) => {
        const rootNodes = [
            'waveforms-root',
            'timings-root',
            'antennas-root',
            'platforms-root',
        ];
        if (nodeId && rootNodes.includes(nodeId)) {
            return;
        }
        selectItem(nodeId);
    };

    return (
        <Box
            sx={{
                height: '100%',
                display: 'flex',
                flexDirection: 'column',
            }}
        >
            {/* 1. Fixed Header */}
            <Box
                sx={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'center',
                    flexShrink: 0,
                    px: 2,
                    pt: 2,
                    pb: 1,
                }}
            >
                <Typography variant="overline" sx={{ color: 'text.secondary' }}>
                    Scenario Explorer
                </Typography>
                <Box>
                    <ScenarioIO />
                </Box>
            </Box>
            <Divider sx={{ mx: 2, mb: 1 }} />

            {/* 2. Scrollable Content Area */}
            <Box
                sx={{
                    flexGrow: 1, // Takes up all remaining space
                    overflowY: 'auto', // Enables vertical scrolling
                    minHeight: 0, // Crucial for flexbox scrolling
                    px: 2, // Apply horizontal padding inside scroll area
                }}
            >
                <SimpleTreeView
                    selectedItems={selectedComponentId ?? selectedItemId}
                    onSelectedItemsChange={handleSelect}
                    slots={{
                        collapseIcon: ExpandMoreIcon,
                        expandIcon: ChevronRightIcon,
                    }}
                    sx={{
                        // Show remove button on hover
                        '& .MuiTreeItem-content .remove-button': {
                            visibility: 'hidden',
                        },
                        '& .MuiTreeItem-content:hover .remove-button': {
                            visibility: 'visible',
                        },
                        '& .MuiTreeItem-content': { py: 0.5 },
                    }}
                >
                    <TreeItem
                        itemId={globalParameters.id}
                        label={
                            <Box sx={{ display: 'flex', alignItems: 'center' }}>
                                <PublicIcon sx={{ mr: 1 }} fontSize="small" />
                                <Typography variant="body2">
                                    Global Parameters
                                </Typography>
                            </Box>
                        }
                    />
                    <TreeItem
                        itemId="waveforms-root"
                        label={
                            <Box
                                sx={{
                                    display: 'flex',
                                    alignItems: 'center',
                                    width: '100%',
                                }}
                            >
                                <WavesIcon sx={{ mr: 1 }} fontSize="small" />
                                <SectionHeader
                                    title="Waveforms"
                                    onAdd={addWaveform}
                                />
                            </Box>
                        }
                    >
                        {waveforms.map((waveform) => (
                            <TreeItem
                                key={waveform.id}
                                itemId={waveform.id}
                                label={
                                    <Box
                                        sx={{
                                            display: 'flex',
                                            alignItems: 'center',
                                            width: '100%',
                                        }}
                                    >
                                        <WavesIcon
                                            sx={{ mr: 1 }}
                                            fontSize="small"
                                        />
                                        <Typography
                                            variant="body2"
                                            sx={{
                                                flexGrow: 1,
                                                whiteSpace: 'nowrap',
                                                overflow: 'hidden',
                                                textOverflow: 'ellipsis',
                                                minWidth: 0,
                                                pr: 1,
                                            }}
                                        >
                                            {waveform.name}
                                        </Typography>
                                        <IconButton
                                            size="small"
                                            className="remove-button"
                                            onClick={(e) => {
                                                e.stopPropagation();
                                                removeItem(waveform.id);
                                            }}
                                        >
                                            <RemoveIcon fontSize="inherit" />
                                        </IconButton>
                                    </Box>
                                }
                            />
                        ))}
                    </TreeItem>
                    <TreeItem
                        itemId="timings-root"
                        label={
                            <Box
                                sx={{
                                    display: 'flex',
                                    alignItems: 'center',
                                    width: '100%',
                                }}
                            >
                                <TimerIcon sx={{ mr: 1 }} fontSize="small" />
                                <SectionHeader
                                    title="Timings"
                                    onAdd={addTiming}
                                />
                            </Box>
                        }
                    >
                        {timings.map((timing) => (
                            <TreeItem
                                key={timing.id}
                                itemId={timing.id}
                                label={
                                    <Box
                                        sx={{
                                            display: 'flex',
                                            alignItems: 'center',
                                            width: '100%',
                                        }}
                                    >
                                        <TimerIcon
                                            sx={{ mr: 1 }}
                                            fontSize="small"
                                        />
                                        <Typography
                                            variant="body2"
                                            sx={{
                                                flexGrow: 1,
                                                whiteSpace: 'nowrap',
                                                overflow: 'hidden',
                                                textOverflow: 'ellipsis',
                                                minWidth: 0,
                                                pr: 1,
                                            }}
                                        >
                                            {timing.name}
                                        </Typography>
                                        <IconButton
                                            size="small"
                                            className="remove-button"
                                            onClick={(e) => {
                                                e.stopPropagation();
                                                removeItem(timing.id);
                                            }}
                                        >
                                            <RemoveIcon fontSize="inherit" />
                                        </IconButton>
                                    </Box>
                                }
                            />
                        ))}
                    </TreeItem>
                    <TreeItem
                        itemId="antennas-root"
                        label={
                            <Box
                                sx={{
                                    display: 'flex',
                                    alignItems: 'center',
                                    width: '100%',
                                }}
                            >
                                <SettingsInputAntennaIcon
                                    sx={{ mr: 1 }}
                                    fontSize="small"
                                />
                                <SectionHeader
                                    title="Antennas"
                                    onAdd={addAntenna}
                                />
                            </Box>
                        }
                    >
                        {antennas.map((antenna) => (
                            <TreeItem
                                key={antenna.id}
                                itemId={antenna.id}
                                label={
                                    <Box
                                        sx={{
                                            display: 'flex',
                                            alignItems: 'center',
                                            width: '100%',
                                        }}
                                    >
                                        <SettingsInputAntennaIcon
                                            sx={{ mr: 1 }}
                                            fontSize="small"
                                        />
                                        <Typography
                                            variant="body2"
                                            sx={{
                                                flexGrow: 1,
                                                whiteSpace: 'nowrap',
                                                overflow: 'hidden',
                                                textOverflow: 'ellipsis',
                                                minWidth: 0,
                                                pr: 1,
                                            }}
                                        >
                                            {antenna.name}
                                        </Typography>
                                        <IconButton
                                            size="small"
                                            className="remove-button"
                                            onClick={(e) => {
                                                e.stopPropagation();
                                                removeItem(antenna.id);
                                            }}
                                        >
                                            <RemoveIcon fontSize="inherit" />
                                        </IconButton>
                                    </Box>
                                }
                            />
                        ))}
                    </TreeItem>
                    <TreeItem
                        itemId="platforms-root"
                        label={
                            <Box
                                sx={{
                                    display: 'flex',
                                    alignItems: 'center',
                                    width: '100%',
                                }}
                            >
                                <FlightIcon sx={{ mr: 1 }} fontSize="small" />
                                <SectionHeader
                                    title="Platforms"
                                    onAdd={addPlatform}
                                />
                            </Box>
                        }
                    >
                        {platforms.map((platform) => {
                            return (
                                <TreeItem
                                    key={platform.id}
                                    itemId={platform.id}
                                    label={
                                        <Box
                                            sx={{
                                                display: 'flex',
                                                alignItems: 'center',
                                                width: '100%',
                                            }}
                                        >
                                            {getPlatformIcon(platform)}
                                            <Typography
                                                variant="body2"
                                                sx={{
                                                    flexGrow: 1,
                                                    whiteSpace: 'nowrap',
                                                    overflow: 'hidden',
                                                    textOverflow: 'ellipsis',
                                                    minWidth: 0,
                                                    pr: 1,
                                                }}
                                            >
                                                {platform.name}
                                            </Typography>
                                            <IconButton
                                                size="small"
                                                className="remove-button"
                                                onClick={(e) => {
                                                    e.stopPropagation();
                                                    removeItem(platform.id);
                                                }}
                                            >
                                                <RemoveIcon fontSize="inherit" />
                                            </IconButton>
                                        </Box>
                                    }
                                >
                                    {platform.components.map((component) => (
                                        <TreeItem
                                            key={component.id}
                                            itemId={component.id}
                                            label={
                                                <Box
                                                    sx={{
                                                        display: 'flex',
                                                        alignItems: 'center',
                                                        width: '100%',
                                                    }}
                                                >
                                                    {getComponentIcon(
                                                        component
                                                    )}
                                                    <Typography
                                                        variant="body2"
                                                        sx={{
                                                            flexGrow: 1,
                                                            whiteSpace:
                                                                'nowrap',
                                                            overflow: 'hidden',
                                                            textOverflow:
                                                                'ellipsis',
                                                            minWidth: 0,
                                                            pr: 1,
                                                        }}
                                                    >
                                                        {component.name}
                                                    </Typography>
                                                    <IconButton
                                                        size="small"
                                                        className="remove-button"
                                                        onClick={(e) => {
                                                            e.stopPropagation();
                                                            removePlatformComponent(
                                                                platform.id,
                                                                component.id
                                                            );
                                                        }}
                                                    >
                                                        <RemoveIcon fontSize="inherit" />
                                                    </IconButton>
                                                </Box>
                                            }
                                        />
                                    ))}
                                </TreeItem>
                            );
                        })}
                    </TreeItem>
                </SimpleTreeView>
            </Box>
        </Box>
    );
}
