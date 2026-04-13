// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import CenterFocusStrongIcon from '@mui/icons-material/CenterFocusStrong';
import CloseIcon from '@mui/icons-material/Close';

import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import TravelExploreIcon from '@mui/icons-material/TravelExplore';
import TuneIcon from '@mui/icons-material/Tune';
import VideocamIcon from '@mui/icons-material/Videocam';
import VideocamOffIcon from '@mui/icons-material/VideocamOff';
import {
    Accordion,
    AccordionDetails,
    AccordionSummary,
    Box,
    Checkbox,
    Collapse,
    Divider,
    FormControlLabel,
    IconButton,
    Paper,
    Stack,
    Switch,
    Tooltip,
    Typography,
} from '@mui/material';
import { useState } from 'react';
import { useScenarioStore, VisualizationLayers } from '@/stores/scenarioStore';

// TODO: users should be able to drag the view controls button or pane to relocate it in the preview area. It should dynamically determine whether to expand upwards or downwards depending on whether it is opening from the bottom of the preview area or the top of the preview area.

export default function ViewControls() {
    const frameScene = useScenarioStore((s) => s.frameScene);
    const focusOnItem = useScenarioStore((s) => s.focusOnItem);
    const toggleFollowItem = useScenarioStore((s) => s.toggleFollowItem);
    const selectedItemId = useScenarioStore((s) => s.selectedItemId);
    const viewControlAction = useScenarioStore((s) => s.viewControlAction);
    const visibility = useScenarioStore((s) => s.visibility);
    const toggleLayer = useScenarioStore((s) => s.toggleLayer);

    const [expanded, setExpanded] = useState(true);

    // Helper to determine if we are following the currently selected item
    const isFollowing =
        viewControlAction.type === 'follow' &&
        viewControlAction.targetId === selectedItemId;

    const handleToggle = (key: keyof VisualizationLayers) => {
        toggleLayer(key);
    };

    if (!expanded) {
        return (
            <Tooltip title="Open View Controls" placement="right">
                <Paper
                    elevation={4}
                    sx={{
                        width: 40,
                        height: 40,
                        display: 'flex',
                        alignItems: 'center',
                        justifyContent: 'center',
                        borderRadius: '50%',
                        overflow: 'hidden',
                    }}
                >
                    <IconButton onClick={() => setExpanded(true)} size="small">
                        <TuneIcon fontSize="small" />
                    </IconButton>
                </Paper>
            </Tooltip>
        );
    }

    return (
        <Paper
            elevation={4}
            sx={{
                width: 280,
                maxHeight: '80vh',
                display: 'flex',
                flexDirection: 'column',
                overflow: 'hidden',
                borderRadius: 2,
                backgroundColor: 'background.paper',
            }}
        >
            {/* 1. Header & Camera Controls (Always Visible) */}
            <Box
                sx={{
                    p: 1,
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'space-between',
                    bgcolor: 'action.hover',
                    borderBottom: 1,
                    borderColor: 'divider',
                }}
            >
                <Stack direction="row" spacing={1}>
                    <Tooltip title="Frame Scene">
                        <IconButton size="small" onClick={frameScene}>
                            <TravelExploreIcon fontSize="small" />
                        </IconButton>
                    </Tooltip>
                    <Tooltip title="Focus on Selected">
                        <span>
                            <IconButton
                                size="small"
                                onClick={() =>
                                    selectedItemId &&
                                    focusOnItem(selectedItemId)
                                }
                                disabled={!selectedItemId}
                                color="primary"
                            >
                                <CenterFocusStrongIcon fontSize="small" />
                            </IconButton>
                        </span>
                    </Tooltip>
                    <Tooltip
                        title={
                            isFollowing ? 'Stop Following' : 'Follow Selected'
                        }
                    >
                        <span>
                            <IconButton
                                size="small"
                                onClick={() =>
                                    selectedItemId &&
                                    toggleFollowItem(selectedItemId)
                                }
                                disabled={!selectedItemId}
                                color={isFollowing ? 'secondary' : 'default'}
                                sx={{
                                    bgcolor: isFollowing
                                        ? 'secondary.main'
                                        : 'transparent',
                                    color: isFollowing ? 'white' : 'inherit',
                                    '&:hover': {
                                        bgcolor: isFollowing
                                            ? 'secondary.dark'
                                            : 'action.hover',
                                    },
                                }}
                            >
                                {isFollowing ? (
                                    <VideocamIcon fontSize="small" />
                                ) : (
                                    <VideocamOffIcon fontSize="small" />
                                )}
                            </IconButton>
                        </span>
                    </Tooltip>
                </Stack>
                <IconButton size="small" onClick={() => setExpanded(false)}>
                    <CloseIcon fontSize="small" />
                </IconButton>
            </Box>

            {/* 2. Scrollable Settings Area */}
            <Box sx={{ overflowY: 'auto', p: 0 }}>
                {/* SECTION: SCENE OBJECTS */}
                <Accordion defaultExpanded disableGutters elevation={0}>
                    <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                        <Typography variant="subtitle2">
                            Scene Objects
                        </Typography>
                    </AccordionSummary>
                    <AccordionDetails sx={{ pt: 0, pb: 2 }}>
                        <Stack spacing={0.5}>
                            <FormControlLabel
                                control={
                                    <Switch
                                        size="small"
                                        checked={visibility.showPlatforms}
                                        onChange={() =>
                                            handleToggle('showPlatforms')
                                        }
                                    />
                                }
                                label={
                                    <Typography variant="body2">
                                        Platforms
                                    </Typography>
                                }
                            />
                            <Collapse in={visibility.showPlatforms}>
                                <Box
                                    sx={{
                                        pl: 4,
                                        display: 'flex',
                                        flexDirection: 'column',
                                    }}
                                >
                                    <FormControlLabel
                                        control={
                                            <Checkbox
                                                size="small"
                                                checked={
                                                    visibility.showPlatformLabels
                                                }
                                                onChange={() =>
                                                    handleToggle(
                                                        'showPlatformLabels'
                                                    )
                                                }
                                            />
                                        }
                                        label={
                                            <Typography variant="caption">
                                                Labels
                                            </Typography>
                                        }
                                    />
                                    <FormControlLabel
                                        control={
                                            <Checkbox
                                                size="small"
                                                checked={visibility.showAxes}
                                                onChange={() =>
                                                    handleToggle('showAxes')
                                                }
                                            />
                                        }
                                        label={
                                            <Typography variant="caption">
                                                Body Axes
                                            </Typography>
                                        }
                                    />
                                    <FormControlLabel
                                        control={
                                            <Checkbox
                                                size="small"
                                                checked={
                                                    visibility.showVelocities
                                                }
                                                onChange={() =>
                                                    handleToggle(
                                                        'showVelocities'
                                                    )
                                                }
                                            />
                                        }
                                        label={
                                            <Typography variant="caption">
                                                Velocity Vectors
                                            </Typography>
                                        }
                                    />
                                    <FormControlLabel
                                        control={
                                            <Checkbox
                                                size="small"
                                                checked={
                                                    visibility.showPatterns
                                                }
                                                onChange={() =>
                                                    handleToggle('showPatterns')
                                                }
                                            />
                                        }
                                        label={
                                            <Typography variant="caption">
                                                Antenna Patterns
                                            </Typography>
                                        }
                                    />
                                    <FormControlLabel
                                        control={
                                            <Checkbox
                                                size="small"
                                                checked={
                                                    visibility.showBoresights
                                                }
                                                onChange={() =>
                                                    handleToggle(
                                                        'showBoresights'
                                                    )
                                                }
                                            />
                                        }
                                        label={
                                            <Typography variant="caption">
                                                Boresight Arrows
                                            </Typography>
                                        }
                                    />
                                </Box>
                            </Collapse>

                            <Divider sx={{ my: 1 }} />

                            <FormControlLabel
                                control={
                                    <Switch
                                        size="small"
                                        checked={visibility.showMotionPaths}
                                        onChange={() =>
                                            handleToggle('showMotionPaths')
                                        }
                                    />
                                }
                                label={
                                    <Typography variant="body2">
                                        Motion Paths
                                    </Typography>
                                }
                            />
                        </Stack>
                    </AccordionDetails>
                </Accordion>

                <Divider />

                {/* SECTION: RF LINKS */}
                <Accordion defaultExpanded disableGutters elevation={0}>
                    <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                        <Typography variant="subtitle2">RF Links</Typography>
                    </AccordionSummary>
                    <AccordionDetails sx={{ pt: 0 }}>
                        <Stack spacing={0.5}>
                            <FormControlLabel
                                control={
                                    <Switch
                                        size="small"
                                        checked={visibility.showLinks}
                                        onChange={() =>
                                            handleToggle('showLinks')
                                        }
                                    />
                                }
                                label={
                                    <Typography variant="body2">
                                        Show Links
                                    </Typography>
                                }
                            />
                            <Collapse in={visibility.showLinks}>
                                <Box sx={{ pl: 1 }}>
                                    <FormControlLabel
                                        control={
                                            <Checkbox
                                                size="small"
                                                checked={
                                                    visibility.showLinkLabels
                                                }
                                                onChange={() =>
                                                    handleToggle(
                                                        'showLinkLabels'
                                                    )
                                                }
                                            />
                                        }
                                        label={
                                            <Typography variant="body2">
                                                Show Labels
                                            </Typography>
                                        }
                                    />
                                    <Typography
                                        variant="caption"
                                        color="text.secondary"
                                        sx={{ mt: 1, display: 'block' }}
                                    >
                                        Filter by Type
                                    </Typography>
                                    <Box
                                        sx={{
                                            pl: 1,
                                            display: 'flex',
                                            flexDirection: 'column',
                                        }}
                                    >
                                        <FormControlLabel
                                            control={
                                                <Checkbox
                                                    size="small"
                                                    checked={
                                                        visibility.showLinkMonostatic
                                                    }
                                                    onChange={() =>
                                                        handleToggle(
                                                            'showLinkMonostatic'
                                                        )
                                                    }
                                                />
                                            }
                                            label={
                                                <Typography variant="caption">
                                                    Monostatic
                                                </Typography>
                                            }
                                        />
                                        <FormControlLabel
                                            control={
                                                <Checkbox
                                                    size="small"
                                                    checked={
                                                        visibility.showLinkIlluminator
                                                    }
                                                    onChange={() =>
                                                        handleToggle(
                                                            'showLinkIlluminator'
                                                        )
                                                    }
                                                />
                                            }
                                            label={
                                                <Typography variant="caption">
                                                    Bistatic (Illuminator)
                                                </Typography>
                                            }
                                        />
                                        <FormControlLabel
                                            control={
                                                <Checkbox
                                                    size="small"
                                                    checked={
                                                        visibility.showLinkScattered
                                                    }
                                                    onChange={() =>
                                                        handleToggle(
                                                            'showLinkScattered'
                                                        )
                                                    }
                                                />
                                            }
                                            label={
                                                <Typography variant="caption">
                                                    Bistatic (Scattered)
                                                </Typography>
                                            }
                                        />
                                        <FormControlLabel
                                            control={
                                                <Checkbox
                                                    size="small"
                                                    checked={
                                                        visibility.showLinkDirect
                                                    }
                                                    onChange={() =>
                                                        handleToggle(
                                                            'showLinkDirect'
                                                        )
                                                    }
                                                />
                                            }
                                            label={
                                                <Typography variant="caption">
                                                    Direct (Interference)
                                                </Typography>
                                            }
                                        />
                                    </Box>
                                </Box>
                            </Collapse>
                        </Stack>
                    </AccordionDetails>
                </Accordion>
            </Box>
        </Paper>
    );
}
