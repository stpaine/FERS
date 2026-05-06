// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import AddIcon from '@mui/icons-material/Add';
import DeleteIcon from '@mui/icons-material/Delete';
import EditIcon from '@mui/icons-material/Edit';
import {
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    FormControl,
    IconButton,
    InputLabel,
    MenuItem,
    Select,
    Typography,
} from '@mui/material';
import { useState } from 'react';
import {
    Platform,
    PlatformComponent,
    PositionWaypoint,
    RotationWaypoint,
    useScenarioStore,
} from '@/stores/scenarioStore';
import { generateSimId } from '@/stores/scenarioStore/idUtils';
import { BufferedTextField, NumberField, Section } from './InspectorControls';
import { PlatformComponentInspector } from './PlatformComponentInspector';

interface PlatformInspectorProps {
    item: Platform;
    selectedComponentId: string | null;
}

type InterpolationType = 'static' | 'linear' | 'cubic';

interface WaypointEditDialogProps {
    open: boolean;
    onClose: (
        finalWaypoint: PositionWaypoint | RotationWaypoint | null
    ) => void;
    waypoint: PositionWaypoint | RotationWaypoint | null;
    waypointType: 'position' | 'rotation';
    angleUnitLabel: 'deg' | 'rad';
}

const createDefaultPositionWaypoint = (): PositionWaypoint => ({
    id: generateSimId('Platform'),
    x: 0,
    y: 0,
    altitude: 0,
    time: 0,
});

const createDefaultRotationWaypoint = (): RotationWaypoint => ({
    id: generateSimId('Platform'),
    azimuth: 0,
    elevation: 0,
    time: 0,
});

export function ensureCubicPositionWaypoints(
    waypoints: PositionWaypoint[]
): PositionWaypoint[] {
    if (waypoints.length >= 2) {
        return waypoints;
    }

    const first = waypoints[0] ?? createDefaultPositionWaypoint();
    return [
        first,
        {
            ...first,
            id: generateSimId('Platform'),
            time: first.time + 1,
        },
    ];
}

export function ensureCubicRotationWaypoints(
    waypoints: RotationWaypoint[]
): RotationWaypoint[] {
    if (waypoints.length >= 2) {
        return waypoints;
    }

    const first = waypoints[0] ?? createDefaultRotationWaypoint();
    return [
        first,
        {
            ...first,
            id: generateSimId('Platform'),
            time: first.time + 1,
        },
    ];
}

function WaypointEditDialog({
    open,
    onClose,
    waypoint,
    waypointType,
    angleUnitLabel,
}: WaypointEditDialogProps) {
    const [editedWaypoint, setEditedWaypoint] = useState(waypoint);

    const handleFieldChange = (field: string, value: number | null) => {
        setEditedWaypoint((prev) => {
            if (!prev) return null;
            return { ...prev, [field]: value };
        });
    };

    const handleClose = () => {
        if (!editedWaypoint) {
            onClose(null);
            return;
        }

        // Create a mutable copy and cast to a record for safe dynamic access.
        const sanitizedWaypoint: Record<string, unknown> = {
            ...editedWaypoint,
        };

        // Sanitize the local state on close: convert any nulls back to 0.
        for (const key in sanitizedWaypoint) {
            if (
                Object.hasOwn(sanitizedWaypoint, key) &&
                sanitizedWaypoint[key] === null
            ) {
                sanitizedWaypoint[key] = 0;
            }
        }
        // Cast back to the original type before calling the parent callback.
        onClose(sanitizedWaypoint as PositionWaypoint | RotationWaypoint);
    };

    if (!editedWaypoint) return null;

    return (
        <Dialog open={open} onClose={handleClose}>
            <DialogTitle>Edit Waypoint</DialogTitle>
            <DialogContent>
                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'column',
                        gap: 2,
                        pt: 1,
                    }}
                >
                    {waypointType === 'position' && 'x' in editedWaypoint && (
                        <>
                            <NumberField
                                label="X"
                                value={editedWaypoint.x}
                                emptyBehavior="revert"
                                onChange={(v) => handleFieldChange('x', v)}
                            />
                            <NumberField
                                label="Y"
                                value={editedWaypoint.y}
                                emptyBehavior="revert"
                                onChange={(v) => handleFieldChange('y', v)}
                            />
                            <NumberField
                                label="Altitude"
                                value={editedWaypoint.altitude}
                                emptyBehavior="revert"
                                onChange={(v) =>
                                    handleFieldChange('altitude', v)
                                }
                            />
                        </>
                    )}
                    {waypointType === 'rotation' &&
                        'azimuth' in editedWaypoint && (
                            <>
                                <NumberField
                                    label={`Azimuth (${angleUnitLabel})`}
                                    value={editedWaypoint.azimuth}
                                    emptyBehavior="revert"
                                    onChange={(v) =>
                                        handleFieldChange('azimuth', v)
                                    }
                                />
                                <NumberField
                                    label={`Elevation (${angleUnitLabel})`}
                                    value={editedWaypoint.elevation}
                                    emptyBehavior="revert"
                                    onChange={(v) =>
                                        handleFieldChange('elevation', v)
                                    }
                                />
                            </>
                        )}
                    <NumberField
                        label="Time (s)"
                        value={editedWaypoint.time}
                        emptyBehavior="revert"
                        onChange={(v) => handleFieldChange('time', v)}
                    />
                </Box>
            </DialogContent>
            <DialogActions>
                <Button onClick={handleClose}>Close</Button>
            </DialogActions>
        </Dialog>
    );
}

export function PlatformInspector({
    item,
    selectedComponentId,
}: PlatformInspectorProps) {
    const {
        updateItem,
        addPositionWaypoint,
        removePositionWaypoint,
        addRotationWaypoint,
        removeRotationWaypoint,
        addPlatformComponent,
        removePlatformComponent,
    } = useScenarioStore.getState();
    const angleUnitLabel = useScenarioStore(
        (state) => state.globalParameters.rotationAngleUnit
    );

    const handleChange = (path: string, value: unknown) =>
        updateItem(item.id, path, value);

    const allowMultiplePosWaypoints =
        item.motionPath.interpolation !== 'static';

    const [editingWaypointInfo, setEditingWaypointInfo] = useState<{
        type: 'position' | 'rotation';
        index: number;
    } | null>(null);

    const [newComponentType, setNewComponentType] =
        useState<PlatformComponent['type']>('monostatic');

    const handleDialogClose = (
        finalWaypoint: PositionWaypoint | RotationWaypoint | null
    ) => {
        if (finalWaypoint && editingWaypointInfo) {
            const { type, index } = editingWaypointInfo;
            const pathPrefix = type === 'position' ? 'motionPath' : 'rotation';
            handleChange(`${pathPrefix}.waypoints.${index}`, finalWaypoint);
        }
        setEditingWaypointInfo(null);
    };

    const currentEditingWaypoint = editingWaypointInfo
        ? editingWaypointInfo.type === 'position'
            ? item.motionPath.waypoints[editingWaypointInfo.index]
            : item.rotation.type === 'path'
              ? item.rotation.waypoints[editingWaypointInfo.index]
              : null
        : null;

    const handleRotationTypeChange = (newType: 'fixed' | 'path') => {
        if (newType === 'fixed' && item.rotation.type !== 'fixed') {
            handleChange('rotation', {
                type: 'fixed',
                startAzimuth: 0,
                startElevation: 0,
                azimuthRate: 0,
                elevationRate: 0,
            });
        } else if (newType === 'path' && item.rotation.type !== 'path') {
            handleChange('rotation', {
                type: 'path',
                interpolation: 'static',
                waypoints: [
                    {
                        id: generateSimId('Platform'),
                        azimuth: 0,
                        elevation: 0,
                        time: 0,
                    },
                ],
            });
        }
    };

    const handleMotionInterpolationChange = (
        interpolation: InterpolationType
    ) => {
        handleChange('motionPath', {
            ...item.motionPath,
            interpolation,
            waypoints:
                interpolation === 'cubic'
                    ? ensureCubicPositionWaypoints(item.motionPath.waypoints)
                    : item.motionPath.waypoints,
        });
    };

    const handleRotationInterpolationChange = (
        interpolation: InterpolationType
    ) => {
        if (item.rotation.type !== 'path') {
            return;
        }

        handleChange('rotation', {
            ...item.rotation,
            interpolation,
            waypoints:
                interpolation === 'cubic'
                    ? ensureCubicRotationWaypoints(item.rotation.waypoints)
                    : item.rotation.waypoints,
        });
    };

    return (
        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
            <BufferedTextField
                label="Name"
                variant="outlined"
                size="small"
                fullWidth
                value={item.name}
                allowEmpty={false}
                onChange={(v) => handleChange('name', v)}
            />

            <Section title="Motion Path">
                <FormControl fullWidth size="small">
                    <InputLabel>Interpolation</InputLabel>
                    <Select
                        label="Interpolation"
                        value={item.motionPath.interpolation}
                        onChange={(e) =>
                            handleMotionInterpolationChange(
                                e.target.value as InterpolationType
                            )
                        }
                    >
                        <MenuItem value="static">Static</MenuItem>
                        <MenuItem value="linear">Linear</MenuItem>
                        <MenuItem value="cubic">Cubic</MenuItem>
                    </Select>
                </FormControl>
                {item.motionPath.waypoints
                    .slice(
                        0,
                        allowMultiplePosWaypoints
                            ? item.motionPath.waypoints.length
                            : 1
                    )
                    .map((wp, i) => (
                        <Box
                            key={wp.id}
                            sx={{
                                display: 'flex',
                                alignItems: 'center',
                                justifyContent: 'space-between',
                                p: 1,
                                border: 1,
                                borderColor: 'divider',
                                borderRadius: 1,
                            }}
                        >
                            <Typography
                                variant="body2"
                                sx={{
                                    flexGrow: 1,
                                    mr: 1,
                                    whiteSpace: 'nowrap',
                                    overflow: 'hidden',
                                    textOverflow: 'ellipsis',
                                }}
                            >
                                {`T: ${wp.time}s, X: ${wp.x}, Y: ${wp.y}, Alt: ${wp.altitude}`}
                            </Typography>
                            <Box sx={{ display: 'flex', flexShrink: 0 }}>
                                <IconButton
                                    size="small"
                                    onClick={() =>
                                        setEditingWaypointInfo({
                                            type: 'position',
                                            index: i,
                                        })
                                    }
                                >
                                    <EditIcon fontSize="small" />
                                </IconButton>
                                {allowMultiplePosWaypoints && (
                                    <IconButton
                                        size="small"
                                        onClick={() =>
                                            removePositionWaypoint(
                                                item.id,
                                                wp.id
                                            )
                                        }
                                        disabled={
                                            item.motionPath.waypoints.length <=
                                            1
                                        }
                                    >
                                        <DeleteIcon fontSize="small" />
                                    </IconButton>
                                )}
                            </Box>
                        </Box>
                    ))}
                {allowMultiplePosWaypoints && (
                    <Button
                        onClick={() => addPositionWaypoint(item.id)}
                        size="small"
                    >
                        Add Waypoint
                    </Button>
                )}
            </Section>

            <Section title="Rotation">
                <FormControl fullWidth size="small">
                    <InputLabel>Rotation Type</InputLabel>
                    <Select
                        label="Rotation Type"
                        value={item.rotation.type}
                        onChange={(e) =>
                            handleRotationTypeChange(
                                e.target.value as 'fixed' | 'path'
                            )
                        }
                    >
                        <MenuItem value="fixed">Fixed Rate</MenuItem>
                        <MenuItem value="path">Waypoint Path</MenuItem>
                    </Select>
                </FormControl>
                {item.rotation.type === 'fixed' && (
                    <>
                        <NumberField
                            label={`Start Azimuth (${angleUnitLabel})`}
                            value={item.rotation.startAzimuth}
                            emptyBehavior="revert"
                            onChange={(v) =>
                                handleChange('rotation.startAzimuth', v)
                            }
                        />
                        <NumberField
                            label={`Start Elevation (${angleUnitLabel})`}
                            value={item.rotation.startElevation}
                            emptyBehavior="revert"
                            onChange={(v) =>
                                handleChange('rotation.startElevation', v)
                            }
                        />
                        <NumberField
                            label={`Azimuth Rate (${angleUnitLabel}/s)`}
                            value={item.rotation.azimuthRate}
                            emptyBehavior="revert"
                            onChange={(v) =>
                                handleChange('rotation.azimuthRate', v)
                            }
                        />
                        <NumberField
                            label={`Elevation Rate (${angleUnitLabel}/s)`}
                            value={item.rotation.elevationRate}
                            emptyBehavior="revert"
                            onChange={(v) =>
                                handleChange('rotation.elevationRate', v)
                            }
                        />
                    </>
                )}
                {item.rotation.type === 'path' &&
                    (() => {
                        const rotation = item.rotation;
                        const allowMultipleWaypoints =
                            rotation.interpolation !== 'static';
                        return (
                            <>
                                <FormControl fullWidth size="small">
                                    <InputLabel>Interpolation</InputLabel>
                                    <Select
                                        label="Interpolation"
                                        value={rotation.interpolation}
                                        onChange={(e) =>
                                            handleRotationInterpolationChange(
                                                e.target
                                                    .value as InterpolationType
                                            )
                                        }
                                    >
                                        <MenuItem value="static">
                                            Static
                                        </MenuItem>
                                        <MenuItem value="linear">
                                            Linear
                                        </MenuItem>
                                        <MenuItem value="cubic">Cubic</MenuItem>
                                    </Select>
                                </FormControl>
                                {rotation.waypoints
                                    .slice(
                                        0,
                                        allowMultipleWaypoints
                                            ? rotation.waypoints.length
                                            : 1
                                    )
                                    .map((wp, i) => (
                                        <Box
                                            key={wp.id}
                                            sx={{
                                                display: 'flex',
                                                alignItems: 'center',
                                                justifyContent: 'space-between',
                                                p: 1,
                                                border: 1,
                                                borderColor: 'divider',
                                                borderRadius: 1,
                                            }}
                                        >
                                            <Typography
                                                variant="body2"
                                                sx={{
                                                    flexGrow: 1,
                                                    mr: 1,
                                                    whiteSpace: 'nowrap',
                                                    overflow: 'hidden',
                                                    textOverflow: 'ellipsis',
                                                }}
                                            >
                                                {`T: ${wp.time}s, Az: ${wp.azimuth} ${angleUnitLabel}, El: ${wp.elevation} ${angleUnitLabel}`}
                                            </Typography>
                                            <Box
                                                sx={{
                                                    display: 'flex',
                                                    flexShrink: 0,
                                                }}
                                            >
                                                <IconButton
                                                    size="small"
                                                    onClick={() =>
                                                        setEditingWaypointInfo({
                                                            type: 'rotation',
                                                            index: i,
                                                        })
                                                    }
                                                >
                                                    <EditIcon fontSize="small" />
                                                </IconButton>
                                                {allowMultipleWaypoints && (
                                                    <IconButton
                                                        size="small"
                                                        onClick={() =>
                                                            removeRotationWaypoint(
                                                                item.id,
                                                                wp.id
                                                            )
                                                        }
                                                        disabled={
                                                            rotation.waypoints
                                                                .length <= 1
                                                        }
                                                    >
                                                        <DeleteIcon fontSize="small" />
                                                    </IconButton>
                                                )}
                                            </Box>
                                        </Box>
                                    ))}
                                {allowMultipleWaypoints && (
                                    <Button
                                        onClick={() =>
                                            addRotationWaypoint(item.id)
                                        }
                                        size="small"
                                    >
                                        Add Waypoint
                                    </Button>
                                )}
                            </>
                        );
                    })()}
            </Section>

            <Section title="Components">
                {item.components.map((comp, index) => (
                    <Box
                        key={comp.id}
                        sx={{
                            border: 1,
                            borderColor: 'divider',
                            borderRadius: 1,
                            p: 2,
                            mb: 2,
                            backgroundColor: 'background.default',
                        }}
                    >
                        <Box
                            sx={{
                                display: 'flex',
                                justifyContent: 'space-between',
                                alignItems: 'center',
                                mb: 2,
                            }}
                        >
                            <Typography variant="subtitle2" color="primary">
                                {comp.type.toUpperCase()}
                                {selectedComponentId === comp.id &&
                                    ' (Selected)'}
                            </Typography>
                            <IconButton
                                size="small"
                                color="error"
                                onClick={() =>
                                    removePlatformComponent(item.id, comp.id)
                                }
                            >
                                <DeleteIcon fontSize="small" />
                            </IconButton>
                        </Box>
                        <PlatformComponentInspector
                            component={comp}
                            platformId={item.id}
                            index={index}
                        />
                    </Box>
                ))}

                <Box
                    sx={{
                        display: 'flex',
                        gap: 1,
                        alignItems: 'center',
                        mt: 1,
                    }}
                >
                    <FormControl size="small" sx={{ flexGrow: 1 }}>
                        <InputLabel>New Component</InputLabel>
                        <Select
                            label="New Component"
                            value={newComponentType}
                            onChange={(e) =>
                                setNewComponentType(
                                    e.target.value as PlatformComponent['type']
                                )
                            }
                        >
                            <MenuItem value="monostatic">
                                Monostatic Radar
                            </MenuItem>
                            <MenuItem value="transmitter">Transmitter</MenuItem>
                            <MenuItem value="receiver">Receiver</MenuItem>
                            <MenuItem value="target">Target</MenuItem>
                        </Select>
                    </FormControl>
                    <Button
                        variant="outlined"
                        onClick={() =>
                            addPlatformComponent(item.id, newComponentType)
                        }
                        startIcon={<AddIcon />}
                    >
                        Add
                    </Button>
                </Box>
            </Section>

            <WaypointEditDialog
                key={
                    editingWaypointInfo
                        ? `${editingWaypointInfo.type}-${editingWaypointInfo.index}`
                        : 'none'
                }
                open={!!editingWaypointInfo}
                onClose={handleDialogClose}
                waypoint={currentEditingWaypoint}
                waypointType={editingWaypointInfo?.type ?? 'position'}
                angleUnitLabel={angleUnitLabel}
            />
        </Box>
    );
}
