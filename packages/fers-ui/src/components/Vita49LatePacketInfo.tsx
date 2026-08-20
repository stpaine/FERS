// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import InfoOutlinedIcon from '@mui/icons-material/InfoOutlined';
import { IconButton, Tooltip } from '@mui/material';

export const latePacketTooltipText = (
    dataPacketCount: number,
    contextPacketCount: number
) =>
    `Packets dispatched more than 1 ms after their scheduled wall-clock deadline. Data: ${dataPacketCount.toLocaleString()}. Context: ${contextPacketCount.toLocaleString()}.`;

export const totalLatePacketCount = (
    dataPacketCount: number,
    contextPacketCount: number
) => dataPacketCount + contextPacketCount;

type Vita49LatePacketInfoProps = {
    dataPacketCount: number;
    contextPacketCount: number;
};

export function Vita49LatePacketInfo({
    dataPacketCount,
    contextPacketCount,
}: Vita49LatePacketInfoProps) {
    const explanation = latePacketTooltipText(
        dataPacketCount,
        contextPacketCount
    );

    return (
        <Tooltip title={explanation} arrow>
            <IconButton
                aria-label={explanation}
                size="small"
                sx={{ ml: 0.25, p: 0.25, color: 'text.secondary' }}
            >
                <InfoOutlinedIcon sx={{ fontSize: 15 }} />
            </IconButton>
        </Tooltip>
    );
}
