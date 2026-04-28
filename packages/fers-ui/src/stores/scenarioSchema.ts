import { z } from 'zod';

// Base numeric type for zod schemas - handles empty strings from forms
const nullableNumber = z.preprocess(
    (val) => (val === '' ? null : val),
    z.number().nullable()
);

// --- SCHEMA DEFINITIONS ---

export const GlobalParametersSchema = z.object({
    id: z.literal('global-parameters'),
    type: z.literal('GlobalParameters'),
    rotationAngleUnit: z.enum(['deg', 'rad']),
    simulation_name: z.string().min(1, 'Simulation name cannot be empty.'),
    start: z.number(),
    end: z.number(),
    rate: z.number().positive('Rate must be positive.'),
    simSamplingRate: nullableNumber.refine((val) => val === null || val > 0, {
        message: 'Sim Sampling Rate must be positive if specified.',
    }),
    c: z.number().positive('Speed of light must be positive.'),
    random_seed: nullableNumber.pipe(z.number().int().nullable()),
    adc_bits: z.number().int().min(0, 'ADC bits cannot be negative.'),
    oversample_ratio: z
        .number()
        .int()
        .min(1, 'Oversample ratio must be at least 1.'),
    origin: z.object({
        latitude: z.number().min(-90).max(90),
        longitude: z.number().min(-180).max(180),
        altitude: z.number(),
    }),
    coordinateSystem: z
        .object({
            frame: z.enum(['ENU', 'UTM', 'ECEF']),
            zone: z.number().int().optional(),
            hemisphere: z.enum(['N', 'S']).optional(),
        })
        .refine(
            (data) => {
                if (data.frame === 'UTM') {
                    return (
                        data.zone !== undefined && data.hemisphere !== undefined
                    );
                }
                return true;
            },
            { message: 'UTM frame requires a zone and hemisphere.' }
        ),
});

const SimIdSchema = z.string().regex(/^\d+$/, 'ID must be a numeric string.');

const BaseWaveformSchema = z.object({
    id: SimIdSchema,
    type: z.literal('Waveform'),
    name: z.string().min(1, 'Waveform name cannot be empty.'),
    power: z.number().min(0, 'Power cannot be negative.'),
    carrier_frequency: z
        .number()
        .positive('Carrier frequency must be positive.'),
});

export const WaveformSchema = z
    .discriminatedUnion('waveformType', [
        BaseWaveformSchema.extend({
            waveformType: z.literal('pulsed_from_file'),
            filename: z
                .string()
                .min(1, 'A filename is required for this waveform type.'),
        }),
        BaseWaveformSchema.extend({
            waveformType: z.literal('cw'),
        }),
        BaseWaveformSchema.extend({
            waveformType: z.literal('fmcw_linear_chirp'),
            direction: z.enum(['up', 'down']),
            chirp_bandwidth: z
                .number()
                .positive('Chirp bandwidth must be positive.'),
            chirp_duration: z
                .number()
                .positive('Chirp duration must be positive.'),
            chirp_period: z.number().positive('Chirp period must be positive.'),
            start_frequency_offset: nullableNumber.pipe(
                z.number().finite().nullable()
            ),
            chirp_count: nullableNumber.pipe(
                z.number().int().positive().nullable()
            ),
        }),
    ])
    .superRefine((data, ctx) => {
        if (
            data.waveformType === 'fmcw_linear_chirp' &&
            data.chirp_period < data.chirp_duration
        ) {
            ctx.addIssue({
                code: 'custom',
                message:
                    'Chirp period must be greater than or equal to chirp duration.',
                path: ['chirp_period'],
            });
        }
    });

export const NoiseEntrySchema = z.object({
    id: SimIdSchema,
    alpha: z.number(),
    weight: z.number(),
});

export const TimingSchema = z.object({
    id: SimIdSchema,
    type: z.literal('Timing'),
    name: z.string().min(1, 'Timing name cannot be empty.'),
    frequency: z.number().positive('Frequency must be positive.'),
    freqOffset: nullableNumber,
    randomFreqOffsetStdev: nullableNumber.pipe(z.number().min(0).nullable()),
    phaseOffset: nullableNumber,
    randomPhaseOffsetStdev: nullableNumber.pipe(z.number().min(0).nullable()),
    noiseEntries: z.array(NoiseEntrySchema),
});

const BaseAntennaSchema = z.object({
    id: SimIdSchema,
    type: z.literal('Antenna'),
    name: z.string().min(1, 'Antenna name cannot be empty.'),
    efficiency: nullableNumber.pipe(z.number().min(0).max(1).nullable()),
    meshScale: nullableNumber.pipe(z.number().positive().nullable()).optional(),
    design_frequency: nullableNumber
        .pipe(z.number().positive().nullable())
        .optional(),
});

export const AntennaSchema = z.discriminatedUnion('pattern', [
    BaseAntennaSchema.extend({ pattern: z.literal('isotropic') }),
    BaseAntennaSchema.extend({
        pattern: z.literal('sinc'),
        alpha: nullableNumber.pipe(z.number().nullable()),
        beta: nullableNumber.pipe(z.number().nullable()),
        gamma: nullableNumber.pipe(z.number().nullable()),
    }),
    BaseAntennaSchema.extend({
        pattern: z.literal('gaussian'),
        azscale: nullableNumber.pipe(z.number().nullable()),
        elscale: nullableNumber.pipe(z.number().nullable()),
    }),
    BaseAntennaSchema.extend({
        pattern: z.literal('squarehorn'),
        diameter: nullableNumber.pipe(z.number().positive().nullable()),
    }),
    BaseAntennaSchema.extend({
        pattern: z.literal('parabolic'),
        diameter: nullableNumber.pipe(z.number().positive().nullable()),
    }),
    BaseAntennaSchema.extend({
        pattern: z.literal('xml'),
        filename: z
            .string()
            .min(1, 'Filename is required for XML pattern.')
            .optional(),
    }),
    BaseAntennaSchema.extend({
        pattern: z.literal('file'),
        filename: z
            .string()
            .min(1, 'Filename is required for file pattern.')
            .optional(),
    }),
]);

export const PositionWaypointSchema = z.object({
    id: SimIdSchema,
    x: z.number(),
    y: z.number(),
    altitude: z.number(),
    time: z.number().min(0, 'Time cannot be negative.'),
});

export const MotionPathSchema = z.object({
    interpolation: z.enum(['static', 'linear', 'cubic']),
    waypoints: z
        .array(PositionWaypointSchema)
        .min(1, 'At least one waypoint is required.'),
});

export const FixedRotationSchema = z.object({
    type: z.literal('fixed'),
    startAzimuth: z.number(),
    startElevation: z.number(),
    azimuthRate: z.number(),
    elevationRate: z.number(),
});

export const RotationWaypointSchema = z.object({
    id: SimIdSchema,
    azimuth: z.number(),
    elevation: z.number(),
    time: z.number().min(0, 'Time cannot be negative.'),
});

export const RotationPathSchema = z.object({
    type: z.literal('path'),
    interpolation: z.enum(['static', 'linear', 'cubic']),
    waypoints: z
        .array(RotationWaypointSchema)
        .min(1, 'At least one waypoint is required.'),
});

export const SchedulePeriodSchema = z.object({
    start: z.number().min(0, 'Start time cannot be negative.'),
    end: z.number().min(0, 'End time cannot be negative.'),
});

const MonostaticComponentSchema = z.object({
    id: SimIdSchema,
    type: z.literal('monostatic'),
    name: z.string().min(1),
    txId: SimIdSchema,
    rxId: SimIdSchema,
    radarType: z.enum(['pulsed', 'cw', 'fmcw']),
    window_skip: nullableNumber,
    window_length: nullableNumber,
    prf: nullableNumber,
    antennaId: SimIdSchema.nullable(),
    waveformId: SimIdSchema.nullable(),
    timingId: SimIdSchema.nullable(),
    noiseTemperature: nullableNumber.pipe(z.number().min(0).nullable()),
    noDirectPaths: z.boolean(),
    noPropagationLoss: z.boolean(),
    schedule: z.array(SchedulePeriodSchema).default([]),
});

const TransmitterComponentSchema = z.object({
    id: SimIdSchema,
    type: z.literal('transmitter'),
    name: z.string().min(1),
    radarType: z.enum(['pulsed', 'cw', 'fmcw']),
    prf: nullableNumber,
    antennaId: SimIdSchema.nullable(),
    waveformId: SimIdSchema.nullable(),
    timingId: SimIdSchema.nullable(),
    schedule: z.array(SchedulePeriodSchema).default([]),
});

const ReceiverComponentSchema = z.object({
    id: SimIdSchema,
    type: z.literal('receiver'),
    name: z.string().min(1),
    radarType: z.enum(['pulsed', 'cw', 'fmcw']),
    window_skip: nullableNumber,
    window_length: nullableNumber,
    prf: nullableNumber,
    antennaId: SimIdSchema.nullable(),
    timingId: SimIdSchema.nullable(),
    noiseTemperature: nullableNumber.pipe(z.number().min(0).nullable()),
    noDirectPaths: z.boolean(),
    noPropagationLoss: z.boolean(),
    schedule: z.array(SchedulePeriodSchema).default([]),
});

const TargetComponentSchema = z.object({
    id: SimIdSchema,
    type: z.literal('target'),
    name: z.string().min(1),
    rcs_type: z.enum(['isotropic', 'file']),
    rcs_value: z.number().optional(),
    rcs_filename: z.string().optional(),
    rcs_model: z.enum(['constant', 'chisquare', 'gamma']),
    rcs_k: z.number().optional(),
});

export const PlatformComponentSchema = z.discriminatedUnion('type', [
    MonostaticComponentSchema,
    TransmitterComponentSchema,
    ReceiverComponentSchema,
    TargetComponentSchema,
]);

export const PlatformSchema = z.object({
    id: SimIdSchema,
    type: z.literal('Platform'),
    name: z.string().min(1, 'Platform name cannot be empty.'),
    motionPath: MotionPathSchema,
    rotation: z.union([FixedRotationSchema, RotationPathSchema]),
    components: z.array(PlatformComponentSchema),
});

export const ScenarioDataSchema = z.object({
    globalParameters: GlobalParametersSchema,
    waveforms: z.array(WaveformSchema),
    timings: z.array(TimingSchema),
    antennas: z.array(AntennaSchema),
    platforms: z.array(PlatformSchema),
});
