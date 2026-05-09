# Simulation Pipelines

This page shows the high-level flow of a FERS run. The diagrams are meant to help users understand what happens to a scenario, why different radar modes produce different output files, and where common settings take effect.

## Whole Run

```mermaid
flowchart TD
    A[Scenario XML] --> B[Validate XML structure]
    B --> C[Load waveforms, antennas, timing, platforms, and targets]
    C --> D[Build simulation world]
    D --> E[Create transmit and receive events]
    E --> F[Run event-driven simulation]
    F --> G[Finalize receiver signals]
    G --> H[Write HDF5 result files]
    D --> I[Optional KML export]
```

Key user settings:

- XML structure comes from [[XML Schema Reference]].
- `<rate>`, `<oversample>`, schedules, platform motion, and radar modes affect event generation and signal rendering.
- Receiver names determine HDF5 output file names.

## UI Scenario Workflow

```mermaid
flowchart TD
    A[Scenario Builder] --> B[Edit global parameters, assets, platforms, and components]
    B --> C[3D viewport previews geometry, pointing, paths, patterns, and RF links]
    C --> D[Timeline previews motion over simulation time]
    B --> E[Export .fersxml]
    B --> F[Simulation Run view]
    F --> G[Run simulation]
    F --> H[Generate KML]
    G --> I[Write HDF5 files]
    G --> J[Show output metadata]
    J --> K[Optional metadata JSON export]
```

The UI preview helps catch setup mistakes before a full run, but the generated HDF5 files are the simulation output. Use external analysis tools for signal processing and plots.

## Scenario Loading

```mermaid
flowchart TD
    A[Main .fersxml file] --> B[Expand include files]
    B --> C[Validate DTD and XSD]
    C --> D[Read global parameters]
    D --> E[Read waveforms]
    E --> F[Read timing sources]
    F --> G[Read antennas]
    G --> H[Read platforms and motion]
    H --> I[Attach transmitters, receivers, monostatics, and targets]
    I --> J[Resolve name references]
    J --> K[Ready to run]
```

What can fail here:

- XML elements are in the wrong order.
- A required parameter is missing.
- A referenced waveform, timing source, or antenna name does not exist.
- A waveform file, antenna file, or RCS file cannot be opened.
- A radar mode does not match its waveform type.

## Event-Driven Simulation

FERS does not simply render every object independently from start to finish. It builds a time-ordered set of events and processes them over the simulation interval.

```mermaid
flowchart TD
    A[Start simulation time] --> B[Schedule initial transmitter and receiver events]
    B --> C{Next event time}
    C --> D[Advance platform positions and rotations]
    D --> E[Process continuous streaming signals up to event time]
    E --> F[Handle event]
    F --> G{More events before end time?}
    G -- yes --> C
    G -- no --> H[Flush remaining streaming and IF data]
    H --> I[Finish output files]
```

Events include:

- Pulsed transmitter pulse starts.
- Pulsed receiver window starts and ends.
- CW/FMCW transmitter streaming starts and ends.
- CW/FMCW receiver streaming starts and ends.

Schedules control when these events are created.

## Pulsed Radar Pipeline

Pulsed radar is organized around transmitted pulses and receiver windows.

```mermaid
flowchart TD
    A[Pulsed waveform file] --> B[Transmitter emits pulse at PRF]
    B --> C[Direct path calculation if enabled]
    B --> D[Target reflection calculation]
    D --> E[Range delay, Doppler, antenna gain, RCS, and propagation loss]
    C --> F[Receiver inbox]
    E --> F
    F --> G[Receiver window opens]
    G --> H[Render I/Q samples for window]
    H --> I[Add timing effects and receiver noise]
    I --> J[Normalize and optional ADC quantization]
    J --> K[Write chunk_N_I and chunk_N_Q datasets]
```

Important pulsed settings:

| Setting | Effect |
| --- | --- |
| Transmitter `<prf>` | Controls how often pulses are emitted. |
| Receiver `<prf>` | Controls how often receive windows are opened. |
| `<window_skip>` | Delays the start of each receive window. |
| `<window_length>` | Controls the duration and sample count of each receive window. |
| `<simSamplingRate>` | Controls geometry interpolation used in pulsed response generation. |
| `<adc_bits>` | Enables final quantization. |

Pulsed output is written as one I/Q dataset pair per receiver window.

## CW Streaming Pipeline

CW output is continuous over the receiver's active intervals.

```mermaid
flowchart TD
    A[CW waveform] --> B[Active transmitter interval]
    B --> C[For each receiver sample]
    C --> D[Update platform geometry]
    D --> E[Calculate direct path if enabled]
    D --> F[Calculate target reflections]
    E --> G[Sum received complex sample]
    F --> G
    G --> H[Apply timing effects and noise]
    H --> I[Normalize and optional ADC quantization]
    I --> J[Write I_data and Q_data]
```

Important CW settings:

| Setting | Effect |
| --- | --- |
| `<rate>` | Output sample rate. |
| Schedules | Active transmit and receive intervals. |
| Antenna patterns and platform rotation | Directional gain over time. |
| `nodirect` | Removes direct transmitter-to-receiver path. |
| `nopropagationloss` | Removes path-loss scaling for debugging. |

## FMCW Pipeline Without Built-In Dechirp

Use `dechirp_mode="none"` when you want FERS to write the received FMCW signal and do the dechirp in your own analysis script.

```mermaid
flowchart TD
    A[FMCW chirp waveform] --> B[Active transmitter interval]
    B --> C[Generate received raw baseband samples]
    C --> D[Apply delays, Doppler, gain, RCS, and path loss]
    D --> E[Sum direct and reflected paths]
    E --> F[Add timing effects and receiver noise]
    F --> G[Normalize and optional ADC quantization]
    G --> H[Write raw I_data and Q_data]
    H --> I[External analysis performs dechirp]
```

Use this mode when:

- You want full control over dechirping.
- You are testing a custom processing chain.
- You want to compare FERS output with another FMCW processor.

## FMCW Pipeline With Built-In Dechirp

Use `dechirp_mode="physical"` or `dechirp_mode="ideal"` when you want FERS to write IF output directly.

```mermaid
flowchart TD
    A[FMCW chirp waveform] --> B[Received raw baseband sample]
    C[Dechirp reference] --> D[Mix received signal with reference]
    B --> D
    D --> E[Low-pass IF filtering if configured]
    E --> F[Resample to if_sample_rate if configured]
    F --> G[Add applicable receiver effects]
    G --> H[Normalize and write I_data and Q_data]
```

Reference choices:

| Reference | Use when |
| --- | --- |
| `attached` | Monostatic radar dechirps using its own waveform. |
| `transmitter` | Bistatic receiver dechirps against a named transmitter. |
| `custom` | Receiver dechirps against a named waveform. |

Dechirp modes:

| Mode | Meaning |
| --- | --- |
| `physical` | Includes timing-source effects in the reference relationship. |
| `ideal` | Uses an idealized reference for cleaner beat-signal analysis. |

## Output Finalization

Finalization is where FERS turns simulated receiver voltage samples into HDF5 datasets.

```mermaid
flowchart TD
    A[Rendered receiver samples] --> B[Add receiver thermal noise]
    B --> C[Apply timing phase effects]
    C --> D[Find fullscale value]
    D --> E{adc_bits > 0?}
    E -- yes --> F[Quantize samples]
    E -- no --> G[Normalize floating-point samples]
    F --> H[Write HDF5 datasets]
    G --> H
    H --> I[Write metadata attributes]
```

Important output behavior:

- Stored samples are normalized.
- Use the HDF5 `fullscale` attribute to reconstruct physical I/Q values.
- Pulsed outputs write chunk datasets.
- CW and FMCW outputs write `I_data` and `Q_data`.
- Metadata attributes describe the receiver mode and sampling settings.
