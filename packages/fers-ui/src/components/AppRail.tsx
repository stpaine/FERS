// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import InfoOutlinedIcon from '@mui/icons-material/InfoOutlined';
import PlayArrowIcon from '@mui/icons-material/PlayArrow';
import SettingsIcon from '@mui/icons-material/Settings';
import TerminalIcon from '@mui/icons-material/Terminal';
import ViewInArIcon from '@mui/icons-material/ViewInAr';
import WidgetsIcon from '@mui/icons-material/Widgets';
import {
    Badge,
    Box,
    List,
    ListItemButton,
    ListItemIcon,
    Tooltip,
} from '@mui/material';
import { useFersLogStore } from '@/stores/fersLogStore';

interface AppRailProps {
    activeView: string;
    onViewChange: (view: string) => void;
    onSettingsClick: () => void;
    onAboutClick: () => void;
}

const views = [
    { id: 'scenario', label: 'Scenario Builder', icon: <ViewInArIcon /> },
    { id: 'assets', label: 'Asset Library', icon: <WidgetsIcon /> },
    { id: 'simulation', label: 'Simulation Run', icon: <PlayArrowIcon /> },
];

export default function AppRail({
    activeView,
    onViewChange,
    onSettingsClick,
    onAboutClick,
}: AppRailProps) {
    const logOpen = useFersLogStore((state) => state.isOpen);
    const toggleLogOpen = useFersLogStore((state) => state.toggleOpen);
    const logCount = useFersLogStore((state) => state.entries.length);

    return (
        <Box
            sx={{
                width: 60,
                height: '100%',
                display: 'flex',
                flexDirection: 'column',
                alignItems: 'center',
                borderRight: 1,
                borderColor: 'divider',
                bgcolor: 'background.paper',
                flexShrink: 0,
                overflow: 'hidden',
            }}
        >
            <List sx={{ flexGrow: 1, overflow: 'auto', width: '100%' }}>
                {views.map((view) => (
                    <Tooltip title={view.label} placement="right" key={view.id}>
                        <ListItemButton
                            selected={activeView === view.id}
                            onClick={() => onViewChange(view.id)}
                            sx={{
                                my: 1,
                                justifyContent: 'center',
                                borderRadius: 2,
                                '&.Mui-selected': {
                                    backgroundColor: 'action.selected',
                                },
                            }}
                        >
                            <ListItemIcon sx={{ minWidth: 0 }}>
                                {view.icon}
                            </ListItemIcon>
                        </ListItemButton>
                    </Tooltip>
                ))}
            </List>
            <Box sx={{ pb: 1 }}>
                <Tooltip title="Raw logs" placement="right">
                    <ListItemButton
                        selected={logOpen}
                        onClick={toggleLogOpen}
                        sx={{
                            my: 1,
                            justifyContent: 'center',
                            '&.Mui-selected': {
                                backgroundColor: 'action.selected',
                            },
                        }}
                    >
                        <ListItemIcon sx={{ minWidth: 0 }}>
                            <Badge
                                badgeContent={logCount}
                                max={999}
                                color="secondary"
                                invisible={logCount === 0}
                            >
                                <TerminalIcon />
                            </Badge>
                        </ListItemIcon>
                    </ListItemButton>
                </Tooltip>
                <Tooltip title="Settings" placement="right">
                    <ListItemButton
                        onClick={onSettingsClick}
                        sx={{ my: 1, justifyContent: 'center' }}
                    >
                        <ListItemIcon sx={{ minWidth: 0 }}>
                            <SettingsIcon />
                        </ListItemIcon>
                    </ListItemButton>
                </Tooltip>
                <Tooltip title="About" placement="right">
                    <ListItemButton
                        onClick={onAboutClick}
                        sx={{ my: 1, justifyContent: 'center' }}
                    >
                        <ListItemIcon sx={{ minWidth: 0 }}>
                            <InfoOutlinedIcon />
                        </ListItemIcon>
                    </ListItemButton>
                </Tooltip>
            </Box>
        </Box>
    );
}
