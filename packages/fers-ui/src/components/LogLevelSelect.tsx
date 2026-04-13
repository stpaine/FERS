// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import {
    FormControl,
    InputLabel,
    MenuItem,
    Select,
    type SelectChangeEvent,
} from '@mui/material';
import type { SxProps, Theme } from '@mui/material/styles';
import { invoke } from '@tauri-apps/api/core';
import { useEffect, useState } from 'react';
import {
    type ConfigurableFersLogLevel,
    type FersLogLevel,
    isConfigurableFersLogLevel,
    LOG_LEVEL_OPTIONS,
    useFersLogStore,
} from '@/stores/fersLogStore';
import { useScenarioStore } from '@/stores/scenarioStore';

interface LogLevelSelectProps {
    id: string;
    label?: string;
    sx?: SxProps<Theme>;
}

const logLevelLabel = (level: ConfigurableFersLogLevel) =>
    level === 'OFF' ? 'OFF (disabled)' : level;

export default function LogLevelSelect({
    id,
    label = 'FERS Log Level',
    sx,
}: LogLevelSelectProps) {
    const logLevel = useFersLogStore((state) => state.logLevel);
    const setLogLevel = useFersLogStore((state) => state.setLogLevel);
    const showError = useScenarioStore((state) => state.showError);
    const [isPending, setIsPending] = useState(false);
    const labelId = `${id}-label`;

    useEffect(() => {
        let active = true;

        invoke<FersLogLevel>('get_log_level')
            .then((level) => {
                if (active && isConfigurableFersLogLevel(level)) {
                    setLogLevel(level);
                }
            })
            .catch((error: unknown) => {
                if (active) {
                    showError(`Failed to get FERS log level: ${String(error)}`);
                }
            });

        return () => {
            active = false;
        };
    }, [setLogLevel, showError]);

    const handleChange = async (
        event: SelectChangeEvent<ConfigurableFersLogLevel>
    ) => {
        const nextLevel = event.target.value as ConfigurableFersLogLevel;
        const previousLevel = useFersLogStore.getState().logLevel;

        setLogLevel(nextLevel);
        setIsPending(true);
        try {
            await invoke<void>('set_log_level', { level: nextLevel });
        } catch (error) {
            setLogLevel(previousLevel);
            showError(`Failed to set FERS log level: ${String(error)}`);
        } finally {
            setIsPending(false);
        }
    };

    return (
        <FormControl size="small" sx={sx}>
            <InputLabel id={labelId}>{label}</InputLabel>
            <Select
                labelId={labelId}
                id={id}
                value={logLevel}
                label={label}
                onChange={handleChange}
                disabled={isPending}
            >
                {LOG_LEVEL_OPTIONS.map((level) => (
                    <MenuItem key={level} value={level}>
                        {logLevelLabel(level)}
                    </MenuItem>
                ))}
            </Select>
        </FormControl>
    );
}
