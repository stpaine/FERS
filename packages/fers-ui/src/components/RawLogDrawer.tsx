// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import CloseIcon from '@mui/icons-material/Close';
import DeleteSweepIcon from '@mui/icons-material/DeleteSweep';
import {
    Box,
    Button,
    Divider,
    IconButton,
    TextField,
    Tooltip,
    Typography,
} from '@mui/material';
import type { ChangeEvent, PointerEvent as ReactPointerEvent } from 'react';
import { useEffect, useRef, useState } from 'react';
import {
    type FersLogEntry,
    type FersLogLevel,
    useFersLogStore,
} from '@/stores/fersLogStore';
import LogLevelSelect from './LogLevelSelect';

interface RawLogDrawerProps {
    leftOffset: number;
    width: number;
    onResize: (nextWidth: number) => void;
    onResizeEnd: (nextWidth: number) => void;
}

const levelColor = (level: FersLogLevel) => {
    switch (level) {
        case 'FATAL':
        case 'ERROR':
            return 'error.main';
        case 'WARNING':
            return 'secondary.main';
        case 'TRACE':
        case 'DEBUG':
            return 'text.secondary';
        case 'OFF':
            return 'text.disabled';
        default:
            return 'text.primary';
    }
};

const RawLogLine = ({ entry }: { entry: FersLogEntry }) => (
    <Box
        component="div"
        sx={{
            color: levelColor(entry.level),
            fontFamily: '"Roboto Mono", Consolas, "Liberation Mono", monospace',
            fontSize: 12,
            lineHeight: 1.55,
            whiteSpace: 'pre-wrap',
            wordBreak: 'break-word',
        }}
    >
        {entry.line}
    </Box>
);

export default function RawLogDrawer({
    leftOffset,
    width,
    onResize,
    onResizeEnd,
}: RawLogDrawerProps) {
    const entries = useFersLogStore((state) => state.entries);
    const droppedCount = useFersLogStore((state) => state.droppedCount);
    const maxLines = useFersLogStore((state) => state.maxLines);
    const clearLogs = useFersLogStore((state) => state.clearLogs);
    const setMaxLines = useFersLogStore((state) => state.setMaxLines);
    const setOpen = useFersLogStore((state) => state.setOpen);
    const [dragPointerId, setDragPointerId] = useState<number | null>(null);
    const drawerRef = useRef<HTMLDivElement | null>(null);
    const scrollRef = useRef<HTMLDivElement | null>(null);
    const shouldAutoScrollRef = useRef(true);

    useEffect(() => {
        if (!shouldAutoScrollRef.current || !scrollRef.current) {
            return;
        }

        scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
    }, [entries.length]);

    const handleScroll = () => {
        const element = scrollRef.current;
        if (!element) {
            return;
        }

        shouldAutoScrollRef.current =
            element.scrollHeight - element.scrollTop - element.clientHeight <
            24;
    };

    const handleMaxLinesChange = (event: ChangeEvent<HTMLInputElement>) => {
        const nextValue = Number(event.target.value);
        if (Number.isFinite(nextValue)) {
            setMaxLines(nextValue);
        }
    };

    const getPointerWidth = (clientX: number) => {
        const drawerBounds = drawerRef.current?.getBoundingClientRect();
        if (!drawerBounds) {
            return width;
        }

        return clientX - drawerBounds.left;
    };

    const handleResizePointerDown = (
        event: ReactPointerEvent<HTMLDivElement>
    ) => {
        if (event.pointerType !== 'touch' && event.button !== 0) {
            return;
        }

        event.preventDefault();
        setDragPointerId(event.pointerId);
    };

    useEffect(() => {
        if (dragPointerId === null) {
            return;
        }

        const previousCursor = document.body.style.cursor;
        const previousUserSelect = document.body.style.userSelect;

        const handlePointerMove = (event: PointerEvent) => {
            if (event.pointerId !== dragPointerId) {
                return;
            }

            onResize(getPointerWidth(event.clientX));
        };

        const handlePointerUp = (event: PointerEvent) => {
            if (event.pointerId !== dragPointerId) {
                return;
            }

            onResizeEnd(getPointerWidth(event.clientX));
            setDragPointerId(null);
        };

        document.body.style.cursor = 'col-resize';
        document.body.style.userSelect = 'none';

        window.addEventListener('pointermove', handlePointerMove);
        window.addEventListener('pointerup', handlePointerUp);
        window.addEventListener('pointercancel', handlePointerUp);

        return () => {
            document.body.style.cursor = previousCursor;
            document.body.style.userSelect = previousUserSelect;
            window.removeEventListener('pointermove', handlePointerMove);
            window.removeEventListener('pointerup', handlePointerUp);
            window.removeEventListener('pointercancel', handlePointerUp);
        };
    }, [dragPointerId, onResize, onResizeEnd, width]);

    return (
        <Box
            ref={drawerRef}
            sx={{
                position: 'fixed',
                left: leftOffset,
                top: 0,
                bottom: 0,
                width,
                minWidth: width,
                maxWidth: width,
                bgcolor: 'background.paper',
                borderRight: 1,
                borderColor: 'divider',
                boxShadow: 8,
                zIndex: (theme) => theme.zIndex.drawer,
                display: 'flex',
                flexDirection: 'column',
                flexShrink: 0,
                overflow: 'hidden',
            }}
        >
            <Box
                role="separator"
                aria-orientation="vertical"
                aria-label="Resize raw log viewer"
                onPointerDown={handleResizePointerDown}
                sx={{
                    position: 'absolute',
                    top: 0,
                    right: -5,
                    bottom: 0,
                    width: 10,
                    zIndex: 2,
                    cursor: 'col-resize',
                    touchAction: 'none',
                    '&::before': {
                        content: '""',
                        position: 'absolute',
                        top: 0,
                        bottom: 0,
                        left: '50%',
                        width: 2,
                        borderRadius: 999,
                        transform: 'translateX(-50%)',
                        backgroundColor:
                            dragPointerId === null
                                ? 'rgba(139, 148, 158, 0.35)'
                                : 'primary.main',
                        transition: 'background-color 0.2s ease',
                    },
                    '&:hover::before': {
                        backgroundColor: 'primary.light',
                    },
                }}
            />
            <Box
                sx={{
                    p: 2,
                    display: 'flex',
                    alignItems: 'center',
                    flexWrap: 'wrap',
                    gap: 1.5,
                }}
            >
                <Box sx={{ minWidth: 0, flexGrow: 1 }}>
                    <Typography variant="h6">Raw Logs</Typography>
                    <Typography variant="caption" color="text.secondary">
                        {entries.length} lines retained, drag edge to resize
                    </Typography>
                </Box>
                <LogLevelSelect
                    id="raw-log-level"
                    label="Log level"
                    sx={{ width: 150 }}
                />
                <TextField
                    label="Max lines"
                    type="number"
                    size="small"
                    value={maxLines}
                    onChange={handleMaxLinesChange}
                    sx={{ width: 120 }}
                    slotProps={{
                        htmlInput: {
                            min: 100,
                            max: 20000,
                            step: 100,
                        },
                    }}
                />
                <Tooltip title="Clear logs">
                    <span>
                        <IconButton
                            aria-label="Clear logs"
                            onClick={clearLogs}
                            disabled={entries.length === 0}
                        >
                            <DeleteSweepIcon />
                        </IconButton>
                    </span>
                </Tooltip>
                <Tooltip title="Close logs">
                    <IconButton
                        aria-label="Close logs"
                        onClick={() => setOpen(false)}
                    >
                        <CloseIcon />
                    </IconButton>
                </Tooltip>
            </Box>
            <Divider />
            {droppedCount > 0 && (
                <Box sx={{ px: 2, py: 1, bgcolor: 'action.hover' }}>
                    <Typography variant="caption" color="text.secondary">
                        {droppedCount} older log lines dropped.
                    </Typography>
                </Box>
            )}
            <Box
                ref={scrollRef}
                onScroll={handleScroll}
                sx={{
                    flexGrow: 1,
                    overflow: 'auto',
                    p: 2,
                    bgcolor: 'background.default',
                }}
            >
                {entries.length === 0 ? (
                    <Box
                        sx={{
                            height: '100%',
                            display: 'flex',
                            alignItems: 'center',
                            justifyContent: 'center',
                            textAlign: 'center',
                            px: 3,
                        }}
                    >
                        <Typography variant="body2" color="text.secondary">
                            No FERS logs yet.
                        </Typography>
                    </Box>
                ) : (
                    entries.map((entry) => (
                        <RawLogLine key={entry.sequence} entry={entry} />
                    ))
                )}
            </Box>
            <Divider />
            <Box sx={{ p: 1.5, display: 'flex', justifyContent: 'flex-end' }}>
                <Button size="small" onClick={() => setOpen(false)}>
                    Close
                </Button>
            </Box>
        </Box>
    );
}
