// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { useEffect, useMemo, useRef } from 'react';
import { Box, Typography, Slider, IconButton } from '@mui/material';
import PlayArrowIcon from '@mui/icons-material/PlayArrow';
import PauseIcon from '@mui/icons-material/Pause';
import FastRewindIcon from '@mui/icons-material/FastRewind';
import FastForwardIcon from '@mui/icons-material/FastForward';
import { useScenarioStore } from '@/stores/scenarioStore';

export default function Timeline() {
    const {
        globalParameters,
        isPlaying,
        currentTime,
        togglePlayPause,
        setCurrentTime,
        targetPlaybackDuration,
    } = useScenarioStore();

    const animationFrameRef = useRef(0);
    const lastTimeRef = useRef(0);

    // Calculate dynamic values based on simulation duration and user settings.
    const { simulationDuration, speedFactor, sliderStep, timeStep } =
        useMemo(() => {
            const duration = Math.max(
                0,
                globalParameters.end - globalParameters.start
            );

            let realPlaybackDuration: number;

            if (targetPlaybackDuration !== null && targetPlaybackDuration > 0) {
                realPlaybackDuration = targetPlaybackDuration;
            } else {
                if (duration > 0 && duration < 5) {
                    realPlaybackDuration = 5;
                } else {
                    realPlaybackDuration = duration;
                }
            }

            const factor =
                realPlaybackDuration > 0 ? duration / realPlaybackDuration : 0;
            const sStep = duration > 0 ? duration / 1000 : 0.01;
            const tStep = duration > 0 ? duration / 100 : 0.01;

            return {
                simulationDuration: duration,
                speedFactor: factor,
                sliderStep: sStep,
                timeStep: tStep,
            };
        }, [
            globalParameters.start,
            globalParameters.end,
            targetPlaybackDuration,
        ]);

    // Effect to stop playback when the end is reached.
    useEffect(() => {
        if (isPlaying && currentTime >= globalParameters.end) {
            togglePlayPause();
        }
    }, [isPlaying, currentTime, globalParameters.end, togglePlayPause]);

    useEffect(() => {
        const animate = (now: number) => {
            const deltaTime = (now - lastTimeRef.current) / 1000;
            lastTimeRef.current = now;

            // Let the store handle clamping. The effect above will stop playback.
            setCurrentTime((prevTime) => prevTime + deltaTime * speedFactor);

            animationFrameRef.current = requestAnimationFrame(animate);
        };

        if (isPlaying && simulationDuration > 0) {
            lastTimeRef.current = performance.now();
            animationFrameRef.current = requestAnimationFrame(animate);
        } else {
            cancelAnimationFrame(animationFrameRef.current);
        }

        return () => cancelAnimationFrame(animationFrameRef.current);
    }, [
        isPlaying,
        setCurrentTime,
        simulationDuration,
        speedFactor,
        globalParameters.end,
    ]);

    const handleSliderChange = (_event: Event, newValue: number | number[]) => {
        setCurrentTime(newValue as number);
    };

    return (
        <Box
            sx={{
                display: 'flex',
                alignItems: 'center',
                gap: 2,
                height: '100%',
                px: 2,
                overflow: 'hidden',
            }}
        >
            <Box sx={{ flexShrink: 0 }}>
                <IconButton
                    size="small"
                    onClick={() => setCurrentTime((t) => t - timeStep)}
                >
                    <FastRewindIcon />
                </IconButton>
                <IconButton onClick={togglePlayPause}>
                    {isPlaying ? <PauseIcon /> : <PlayArrowIcon />}
                </IconButton>
                <IconButton
                    size="small"
                    onClick={() => setCurrentTime((t) => t + timeStep)}
                >
                    <FastForwardIcon />
                </IconButton>
            </Box>
            <Typography
                variant="caption"
                sx={{ flexShrink: 0, width: '4em', textAlign: 'center' }}
            >
                {currentTime.toFixed(2)}s
            </Typography>
            <Slider
                value={currentTime}
                onChange={handleSliderChange}
                step={sliderStep}
                min={globalParameters.start}
                max={globalParameters.end}
                sx={{
                    mx: 1,
                    flexGrow: 1,
                    minWidth: 50,
                    // Disable transitions to ensure precise, immediate tracking
                    '& .MuiSlider-thumb, & .MuiSlider-track': {
                        transition: 'none',
                    },
                }}
                valueLabelDisplay="auto"
                valueLabelFormat={(value) => `${value.toFixed(2)}s`}
            />
            <Typography
                variant="caption"
                sx={{ flexShrink: 0, width: '4em', textAlign: 'right' }}
            >
                {globalParameters.end.toFixed(2)}s
            </Typography>
        </Box>
    );
}
