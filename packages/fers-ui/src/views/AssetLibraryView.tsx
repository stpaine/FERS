// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { Box, Typography } from '@mui/material';

export function AssetLibraryView() {
    return (
        <Box sx={{ p: 3 }}>
            <Typography variant="h4">Asset Library</Typography>
            <Typography>
                Manage reusable simulation assets like pulses, antennas, and
                timing sources here.
            </Typography>
        </Box>
    );
}
