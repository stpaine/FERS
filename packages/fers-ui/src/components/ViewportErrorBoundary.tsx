// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import React from 'react';

interface ViewportErrorBoundaryProps {
    children: React.ReactNode;
    fallback?: React.ReactNode;
    onError?: (error: Error) => void;
    resetKey?: string | number;
}

interface ViewportErrorBoundaryState {
    hasError: boolean;
}

export class ViewportErrorBoundary extends React.Component<
    ViewportErrorBoundaryProps,
    ViewportErrorBoundaryState
> {
    override state: ViewportErrorBoundaryState = {
        hasError: false,
    };

    static getDerivedStateFromError(): ViewportErrorBoundaryState {
        return {
            hasError: true,
        };
    }

    override componentDidCatch(error: Error): void {
        this.props.onError?.(error);
    }

    override componentDidUpdate(
        prevProps: Readonly<ViewportErrorBoundaryProps>
    ): void {
        if (this.state.hasError && prevProps.resetKey !== this.props.resetKey) {
            this.setState({ hasError: false });
        }
    }

    override render(): React.ReactNode {
        if (this.state.hasError) {
            return this.props.fallback ?? null;
        }

        return this.props.children;
    }
}
