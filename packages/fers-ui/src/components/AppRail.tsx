// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import {
    Box,
    List,
    ListItemButton,
    ListItemIcon,
    Tooltip,
} from '@mui/material';
import ViewInArIcon from '@mui/icons-material/ViewInAr';
import WidgetsIcon from '@mui/icons-material/Widgets';
import PlayArrowIcon from '@mui/icons-material/PlayArrow';
import BarChartIcon from '@mui/icons-material/BarChart';
import SettingsIcon from '@mui/icons-material/Settings';

interface AppRailProps {
    activeView: string;
    onViewChange: (view: string) => void;
    onSettingsClick: () => void;
}

const views = [
    { id: 'scenario', label: 'Scenario Builder', icon: <ViewInArIcon /> },
    { id: 'assets', label: 'Asset Library', icon: <WidgetsIcon /> },
    { id: 'simulation', label: 'Simulation Run', icon: <PlayArrowIcon /> },
    { id: 'results', label: 'Results Analysis', icon: <BarChartIcon /> },
];

export default function AppRail({
    activeView,
    onViewChange,
    onSettingsClick,
}: AppRailProps) {
    return (
        <Box
            sx={{
                width: 60,
                height: '100%', // Use percentage instead of vh
                display: 'flex',
                flexDirection: 'column',
                alignItems: 'center',
                borderRight: 1,
                borderColor: 'divider',
                bgcolor: 'background.paper',
                flexShrink: 0,
                overflow: 'hidden', // Prevent overflow
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
            </Box>
        </Box>
    );
}
