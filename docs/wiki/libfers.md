# libfers

Most users should run FERS with `fers-cli` or the desktop UI. Use `libfers` only when you want to embed FERS inside another application, test harness, service, or custom user interface.

`libfers` provides the same simulator used by the CLI. It can load scenarios, run simulations, generate KML, report progress, return output metadata, and serialize scenarios.

## When To Use libfers

Use `libfers` when you need to:

- Run FERS from another program instead of a shell command.
- Build your own scenario editor or workflow tool.
- Generate or modify scenarios programmatically.
- Run simulations from an automated test or batch system.
- Receive progress updates while a simulation is running.
- Generate KML or preview data from a loaded scenario.

Do not use `libfers` just to run ordinary `.fersxml` files. The CLI is simpler for that.

## Basic Flow

A program using `libfers` normally does this:

1. Create a FERS context.
2. Configure logging.
3. Load a scenario from a `.fersxml` file or XML string.
4. Choose an output directory.
5. Run the simulation or generate KML.
6. Read the returned metadata if needed.
7. Destroy the context.

The context represents one loaded scenario and its simulation state. Create a separate context when you need an independent scenario.

## Scenario Input

`libfers` can load scenarios from:

- A `.fersxml` file.
- An XML string already held in memory.
- JSON used by FERS UI-style workflows.

File loading is the normal choice. It resolves relative asset paths and processes `<include>` files relative to the scenario file.

XML-string loading is useful when another application already has the XML in memory. It does not process `<include>` files, so embedded applications should expand includes themselves or load from a file.

## Running Simulations

A `libfers` simulation run is synchronous: the call does not return until the simulation is finished or has failed.

If your application has a UI, run the simulation on a background thread and use the progress callback to update the interface.

The simulation writes the same HDF5 files as `fers-cli`. See [[Using FERS]] for output file structure.

HDF5 is the default output backend. Applications can select the FERS VITA 49.2 UDP profile backend at runtime:

```c
fers_enable_vita49_udp_output(context, "127.0.0.1", 4991);
fers_set_vita49_fullscale(context, 1.0);
fers_set_vita49_epoch_unix_nanoseconds(context, 1700000000123456789ULL);
```

The VITA epoch setter is optional. The full-scale setter is required before running in VITA mode. Additional tuning setters are available for maximum UDP payload and sender queue depth. When the queue is full, VITA output blocks the producer until the paced sender frees a slot:

```c
fers_set_vita49_max_udp_payload(context, 1400);
fers_set_vita49_queue_depth(context, 1024);
fers_set_vita49_packet_trace_enabled(context, 1);
```

For the full CLI/API contract, packet layout, metadata fields, and counter semantics, see [[VITA49 Streaming Implementation]].

Do not put output transport selection in `.fersxml`; it remains a scenario/physics format.

## Progress Reporting

Applications can provide a progress callback. Use it for:

- Progress bars.
- Status text.
- Cancel/stop decisions in the calling application.
- Logging in batch systems.

The callback reports high-level simulation progress. It is not intended to report every low-level simulation step.

## Logging And Errors

`libfers` reports logs through the FERS logging system. Applications can:

- Choose a log level.
- Write logs to a file.
- Register a callback to receive formatted log messages.

Most operations return an integer success/failure code. When a call fails, ask `libfers` for the last error message and display it to the user.

Warnings can also be collected after loading or updating a scenario. This is useful for UI workflows where a scenario can be technically loadable but still suspicious.

## Output Metadata

After a run, applications can request JSON metadata describing the output files. This is useful for opening result files without guessing:

- Which files were written.
- Receiver mode.
- Sample rates.
- Sample counts.
- Time span.
- FMCW dechirp and IF settings, when applicable.
- Fullscale information.

Prefer this metadata over hard-coded output assumptions.

When VITA 49.2 output is selected, the metadata includes a `vita49` section with endpoint, epoch, class ID, fixed full-scale, packet sizing, queue depth, receiver mode, and per-stream packet/sample counters when the streaming backend has populated them. Dropped packet/sample counters represent socket send failures, not queue overflow.

## KML Generation

`libfers` can generate a KML file from a loaded scenario. Use this when an application needs a "preview scene on a map" action without running the full simulation.

The KML output depends on the scenario's KML/geospatial settings. Check `<origin>` and `<coordinatesystem>` if exported geometry looks wrong. These settings do not affect the signal simulation.

## Scenario Preview Features

`libfers` includes utility functions for UI-style previews:

- Sample an antenna pattern for plotting.
- Interpolate a platform motion path.
- Interpolate a platform rotation path.
- Calculate preview links between transmitters, receivers, monostatic radars, and targets at a selected time.

These are convenience features for applications. They do not replace a full simulation run.

## Memory Management

Strings and preview data returned by `libfers` are allocated by the library. The calling application is responsible for freeing them with the matching `libfers` free function.

Never use returned pointers after freeing them. Never use a context after destroying it.

## Embedding Checklist

When embedding FERS, make sure your application exposes these user-facing choices:

- Scenario file or generated XML input.
- XML validation on or off. Validation should be on by default.
- Output directory.
- Optional VITA 49.2 UDP endpoint and fixed ADC full-scale.
- Log level and optional log file.
- Worker thread count for larger runs.
- Progress display for long simulations.
- Clear error messages when a scenario cannot be loaded or run.
- A way to open or export the generated HDF5 result files.
- A way to show or save output metadata.

For most applications, the safest design is to treat `.fersxml` as the user-editable scenario format and HDF5 as the analysis output format. Use JSON updates and preview helpers only when building an interactive editor.
