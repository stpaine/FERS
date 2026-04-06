// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

import { z } from 'zod';
import {
    AntennaSchema,
    FixedRotationSchema,
    GlobalParametersSchema,
    MotionPathSchema,
    NoiseEntrySchema,
    PlatformComponentSchema,
    PlatformSchema,
    PositionWaypointSchema,
    RotationPathSchema,
    RotationWaypointSchema,
    ScenarioDataSchema,
    SchedulePeriodSchema,
    TimingSchema,
    WaveformSchema,
} from '../scenarioSchema';

// --- Zod Inferred Types ---
export type GlobalParameters = z.infer<typeof GlobalParametersSchema>;
export type Waveform = z.infer<typeof WaveformSchema>;
export type NoiseEntry = z.infer<typeof NoiseEntrySchema>;
export type Timing = z.infer<typeof TimingSchema>;
export type Antenna = z.infer<typeof AntennaSchema>;
export type PositionWaypoint = z.infer<typeof PositionWaypointSchema>;
export type MotionPath = z.infer<typeof MotionPathSchema>;
export type FixedRotation = z.infer<typeof FixedRotationSchema>;
export type RotationWaypoint = z.infer<typeof RotationWaypointSchema>;
export type RotationPath = z.infer<typeof RotationPathSchema>;
export type PlatformComponent = z.infer<typeof PlatformComponentSchema>;
export type SchedulePeriod = z.infer<typeof SchedulePeriodSchema>;
export type Platform = z.infer<typeof PlatformSchema> & {
    pathPoints?: {
        x: number;
        y: number;
        z: number;
        vx: number;
        vy: number;
        vz: number;
    }[];
    rotationPathPoints?: { azimuth: number; elevation: number }[];
};

// --- Visibility Type Definitions ---
export type VisualizationLayers = {
    showAxes: boolean;
    showPatterns: boolean;
    showBoresights: boolean;
    showLinks: boolean;
    showLinkLabels: boolean;
    showLinkMonostatic: boolean;
    showLinkIlluminator: boolean;
    showLinkScattered: boolean;
    showLinkDirect: boolean;
    showVelocities: boolean;
    showPlatforms: boolean;
    showPlatformLabels: boolean;
    showMotionPaths: boolean;
};

// --- Component Type Definitions ---
export type MonostaticComponent = Extract<
    PlatformComponent,
    { type: 'monostatic' }
>;
export type TransmitterComponent = Extract<
    PlatformComponent,
    { type: 'transmitter' }
>;
export type ReceiverComponent = Extract<
    PlatformComponent,
    { type: 'receiver' }
>;
export type TargetComponent = Extract<PlatformComponent, { type: 'target' }>;

export type ScenarioData = Omit<
    z.infer<typeof ScenarioDataSchema>,
    'platforms'
> & {
    platforms: Platform[];
};

export type ScenarioItem =
    | GlobalParameters
    | Waveform
    | Timing
    | Antenna
    | Platform;

// --- Store Shape ---
export type ViewControlAction = {
    type: 'frame' | 'focus' | 'follow' | null;
    targetId?: string;
    timestamp: number;
};

export type ScenarioState = ScenarioData & {
    selectedItemId: string | null;
    selectedComponentId: string | null;
    isDirty: boolean;
    isPlaying: boolean;
    currentTime: number;
    targetPlaybackDuration: number | null;
    isSimulating: boolean;
    isGeneratingKml: boolean;
    isBackendSyncing: boolean;
    backendVersion: number;
    scenarioFilePath: string | null;
    outputDirectory: string | null;
    antennaPreviewErrors: Record<string, string>;
    errorSnackbar: {
        open: boolean;
        message: string;
    };
    viewControlAction: ViewControlAction;
    visibility: VisualizationLayers;
};

// --- Action Slice Types ---
export type AssetActions = {
    addWaveform: () => void;
    addTiming: () => void;
    addAntenna: () => void;
    addNoiseEntry: (timingId: string) => void;
    removeNoiseEntry: (timingId: string, entryId: string) => void;
    setAntennaPattern: (
        antennaId: string,
        newPattern: Antenna['pattern']
    ) => void;
};

export type PlatformActions = {
    addPlatform: () => void;
    addPositionWaypoint: (platformId: string) => void;
    removePositionWaypoint: (platformId: string, waypointId: string) => void;
    addRotationWaypoint: (platformId: string) => void;
    removeRotationWaypoint: (platformId: string, waypointId: string) => void;
    addPlatformComponent: (
        platformId: string,
        componentType: PlatformComponent['type']
    ) => void;
    removePlatformComponent: (platformId: string, componentId: string) => void;
    setPlatformRcsModel: (
        platformId: string,
        componentId: string,
        newModel: TargetComponent['rcs_model']
    ) => void;
    fetchPlatformPath: (platformId: string) => Promise<void>;
};

export type ScenarioActions = {
    selectItem: (itemId: string | null) => void;
    updateItem: (itemId: string, propertyPath: string, value: unknown) => void;
    setRotationAngleUnit: (
        unit: GlobalParameters['rotationAngleUnit'],
        convertExisting: boolean
    ) => void;
    removeItem: (itemId: string) => void;
    loadScenario: (backendData: unknown) => void;
    resetScenario: () => void;
    setScenarioFilePath: (path: string | null) => void;
    setOutputDirectory: (dir: string | null) => void;
};

export type BackendActions = {
    syncBackend: () => Promise<void>;
    fetchFromBackend: () => Promise<void>;
};

export type PlaybackActions = {
    togglePlayPause: () => void;
    setCurrentTime: (time: number | ((prevTime: number) => number)) => void;
    setTargetPlaybackDuration: (duration: number | null) => void;
    setIsSimulating: (isSimulating: boolean) => void;
    setIsGeneratingKml: (isGeneratingKml: boolean) => void;
};

export type ErrorActions = {
    showError: (message: string) => void;
    hideError: () => void;
    setAntennaPreviewError: (antennaId: string, message: string) => void;
    clearAntennaPreviewError: (antennaId: string) => void;
};

export type ViewControlActions = {
    frameScene: () => void;
    focusOnItem: (itemId: string) => void;
    toggleFollowItem: (itemId: string) => void;
    clearViewControlAction: () => void;
    toggleLayer: (layer: keyof VisualizationLayers) => void;
};

// --- Full Store Type ---
export type FullScenarioActions = AssetActions &
    PlatformActions &
    ScenarioActions &
    BackendActions &
    PlaybackActions &
    ErrorActions &
    ViewControlActions;

export type ScenarioStore = ScenarioState & FullScenarioActions;
