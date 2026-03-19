// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { Box, Typography } from '@mui/material';

export function ResultsView() {
    return (
        <Box sx={{ p: 3 }}>
            <Typography variant="h4">Results Analysis</Typography>
            <Typography>
                Visualize and analyze FERS simulation output data.
            </Typography>
        </Box>
    );
}
