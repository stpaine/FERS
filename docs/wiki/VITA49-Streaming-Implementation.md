# FERS VITA 49.2 UDP Streaming

This page documents the current FERS VITA 49.2 UDP output backend from the source code. It intentionally avoids protocol claims that are not implemented in the repository.

The relevant implementation lives under `packages/libfers/src/serial/vita49/`. CLI option parsing lives in `packages/fers-cli/src/arg_parser.cpp`. Public C API controls live in `packages/libfers/include/libfers/api.h` and `packages/libfers/src/api.cpp`.

## Status

VITA 49.2 UDP output is available, but the project README still marks it as Alpha and Unverified. Treat this page as an implementation contract for the current branch, not as a claim of external interoperability certification.

HDF5 remains the default receiver output backend. VITA output is selected only at runtime through CLI switches or the C API. `.fersxml` scenario files do not select network transport.

## CLI Usage

Minimal VITA UDP run:

```bash
./build/release/packages/fers-cli/fers-cli scenario.fersxml \
  --out-dir=./results \
  --vita49 127.0.0.1:4991 \
  --vita49-fullscale 1.0
```

Deterministic replay epoch:

```bash
./build/release/packages/fers-cli/fers-cli scenario.fersxml \
  --vita49 127.0.0.1:4991 \
  --vita49-fullscale 1.0 \
  --vita49-epoch 1700000000123456789
```

Packet sizing and queue tuning:

```bash
./build/release/packages/fers-cli/fers-cli scenario.fersxml \
  --vita49 127.0.0.1:4991 \
  --vita49-fullscale 1.0 \
  --vita49-max-udp-payload 1400 \
  --vita49-queue-depth 1024
```

CLI VITA options:

| Option | Required | Source-backed behavior |
| --- | --- | --- |
| `--vita49 host:port` | Yes, for VITA mode | Selects VITA UDP output and stores destination host and UDP port. The CLI parser expects exactly `host:port`, rejects an empty host, rejects port `0`, rejects ports above `65535`, and rejects host strings containing another colon. |
| `--vita49-fullscale <positive-real>` | Yes, with `--vita49` | Sets fixed ADC full-scale used to scale complex samples to signed int16 IQ. Must be positive and finite. |
| `--vita49-epoch <unix-nanoseconds>` | No | Sets deterministic VRT UTC epoch as Unix nanoseconds. Must be an unsigned decimal integer that fits the VRT 32-bit UTC seconds field. Maximum accepted value is `4294967295999999999`. |
| `--vita49-max-udp-payload <bytes>` | No | Caps one UDP datagram payload. Valid range is `64..65507` bytes. Default is `1400`. |
| `--vita49-queue-depth <packets>` | No | Sets bounded sender queue depth. Valid range is `1..4294967295`. Default is `1024`. |

When `--vita49` is absent, `fers-cli` uses the default HDF5 backend. When `--vita49` is present, receiver samples are streamed to UDP instead of written through the HDF5 output sink.

## Example Capture Scripts

The repository includes two non-UI VITA examples:

- `examples/vita49_mixed_modes`: pulsed, CW, and FMCW receiver streams in one run.
- `examples/vita49_cw_streaming`: CW streaming capture and analysis.

The capture helpers bind a local UDP socket, run `fers-cli`, and write received datagrams to JSONL with base64 packet payloads.

```bash
export PATH="$PWD/build/release/packages/fers-cli:$PATH"
cd examples/vita49_mixed_modes
python3 capture.py --host 127.0.0.1 --port 4991 --fullscale 1.0
python3 analysis.py
```

`examples/vita49_common.py` decodes the current packet layout. It expects class ID `0xFA52530001000101`, reads packet type from the VRT header, reads signal IQ as big-endian int16 pairs, and reads context ASCII JSON from byte offset `92`.

## C API Usage

Applications using `libfers` can select VITA output before running a simulation:

```c
fers_enable_vita49_udp_output(context, "127.0.0.1", 4991);
fers_set_vita49_fullscale(context, 1.0);
fers_set_vita49_epoch_unix_nanoseconds(context, 1700000000123456789ULL);
fers_set_vita49_max_udp_payload(context, 1400);
fers_set_vita49_queue_depth(context, 1024);
fers_set_vita49_packet_trace_enabled(context, 1);
```

API controls:

| Function | Behavior |
| --- | --- |
| `fers_enable_vita49_udp_output` | Switches the context output mode to VITA UDP and stores destination host/port. Host must be non-null and non-empty. Port must be `1..65535`. |
| `fers_set_vita49_fullscale` | Stores fixed ADC full-scale. Value must be positive and finite. VITA runs are rejected if this is not set to a valid value. |
| `fers_set_vita49_epoch_unix_nanoseconds` | Stores deterministic Unix-nanosecond epoch. Value must fit VRT 32-bit UTC seconds. |
| `fers_set_vita49_max_udp_payload` | Stores maximum UDP payload. Value must be `64..65507`. |
| `fers_set_vita49_queue_depth` | Stores sender queue depth. Value must be greater than zero. |
| `fers_set_vita49_packet_trace_enabled` | Enables or disables packet trace telemetry. Stream counter telemetry is unaffected. This has a C API setter; there is no CLI switch for it. |

Before `fers_run_simulation` or `fers_run_simulation_ex` starts, `libfers` validates the active VITA config. Invalid endpoint, full-scale, payload size, queue depth, or epoch returns an error instead of starting the run.

## Runtime Model

VITA output is implemented at the `core::ReceiverOutputSink` boundary. `core::sim_threading.cpp` chooses `serial::vita49::makeVita49OutputSink(...)` when `core::OutputConfig.mode` is `Vita49Udp`; otherwise it chooses the HDF5 output sink.

The VITA sink creates:

- `Vita49Packetizer`: converts `core::ReceiverSampleBlock` values into serialized VRT packets.
- `PacedSender`: bounded queue and wall-clock pacing thread.
- `UdpSender`: UDP socket sender using `getaddrinfo`, `SOCK_DGRAM`, and `sendto`.

Pacing maps simulation sample time to `std::chrono::steady_clock` time. The sender starts at `params::startTime()`. A packet due time is:

```text
steady_epoch + (packet.first_sample_time - simulation_start_time)
```

If the queue is full, enqueue blocks until the sender frees space. Queue pressure does not create dropped packets. Dropped packet counters are updated when socket send fails.

`finalize()` flushes queued packets, emits close context packets for streams that are still open, stops the sender, and returns final stream stats. `core::sim_threading.cpp` reports `Waiting for VITA output stream drain...` before finalization in VITA mode.

## Stream IDs

FERS allocates one VRT Stream ID per receiver stream descriptor key:

```text
receiver_id + receiver_name + receiver_mode
```

`serial/vita49/stream_registry.cpp` hashes that key with FNV-1a, masks it to a nonzero 31-bit value, and increments on collision until it finds an unused ID. The same key gets the same ID within a run. Do not assume the VRT Stream ID equals the receiver ID.

Tests verify that:

- The same descriptor returns a stable ID.
- Two receivers with colliding low receiver-ID bits do not truncate to the same stream ID.
- The same receiver in different modes, such as pulsed and CW, gets distinct stream IDs.

## Packet Summary

All integer and floating-point fields written by the serializer are big-endian. The serializer writes IEEE 754 binary64 values by bit-casting `RealType` and writing the resulting 64-bit value big-endian.

Current internal class ID:

```text
0xFA52530001000101
```

This value is built from:

- Internal OUI placeholder: `0xFA5253`
- Information class: `0x0001`
- Packet class: `0x0001`
- Version/profile byte: `0x01`

The source comment says this is an internal placeholder until FERS obtains or finalizes an assigned VRT profile identifier.

### Header Fields

`makeHeader(...)` sets these fields:

| Header field | Current source behavior |
| --- | --- |
| Packet type | Signal data uses `0x1`; context uses `0x4`. |
| Class ID present | Set for signal data and context packets. |
| Trailer present | Set for signal data packets; not set for context packets. |
| VITA 49.2 format bit | Set. |
| TSI | UTC (`1`). |
| TSF | Real-time picoseconds (`2`). |
| Packet count | 4-bit counter, rolls over modulo 16. |
| Packet size | 32-bit word count of the whole packet/datagram. |

## Signal Data Packets

Signal data packets contain complex Cartesian IQ samples as interleaved signed 16-bit values:

```text
I0, Q0, I1, Q1, ...
```

Packet byte layout:

| Offset | Size | Field |
| --- | ---: | --- |
| `0` | 4 | VRT header |
| `4` | 4 | Stream ID |
| `8` | 8 | Class ID |
| `16` | 4 | UTC integer seconds |
| `20` | 8 | fractional picoseconds |
| `28` | variable | interleaved big-endian int16 IQ payload |
| last 4 bytes | 4 | trailer |

Fixed overhead is `32` bytes. With the default maximum UDP payload of `1400` bytes, the packetizer can place `342` complex samples in one signal data datagram:

```text
(1400 - 32) / 4 = 342
```

Timestamp behavior:

- The timestamp is the first sample time in the packet.
- `--vita49-epoch` or the configured API epoch supplies the Unix-nanosecond base.
- If no epoch is configured, the sink uses `std::chrono::system_clock::now()` during `initializeRun`.
- Fractional seconds are rounded to picoseconds.
- Integer UTC seconds must fit in a 32-bit unsigned field.

Scaling behavior:

- `--vita49-fullscale` or `fers_set_vita49_fullscale` supplies fixed full-scale.
- Each real and imaginary component is scaled as `round(component / fullscale * 32767)`.
- Exact `+fullscale` maps to `32767`.
- Exact `-fullscale` maps to `-32768`.
- Components greater than `+fullscale` clip to `32767`.
- Components less than `-fullscale` clip to `-32768`.
- `over_range_count` increments once per complex sample where at least one component exceeds the configured full-scale.
- The signal trailer `OverRange` indicator is set when a packet contains at least one clipped complex sample.

Signal trailer indicators are enabled for:

- Calibrated time
- Valid data
- Reference lock
- Over-range
- Sample loss

The corresponding indicator bits are set from the `core::ReceiverSampleBlock` flags, clipping result, and pending sample-loss state.

## Context Packets

The sink emits context packets when `openStream`, `emitContextHeartbeat`, `closeStream`, or finalization close paths call `emitContext`. The current source does not implement a generic "metadata changed" context trigger.

Heartbeat behavior:

- `SimulationEngine` schedules output heartbeats on the simulation clock starting at `params::startTime() + 1.0`.
- `Vita49OutputSink::emitContextHeartbeat` emits a context packet for each open stream when at least one second has elapsed since that stream's last context time.
- The pulsed finalizer path opens a stream, submits rendered sample blocks, and closes the stream; it does not contain its own heartbeat loop.

Context packet byte layout:

| Offset | Size | Field |
| --- | ---: | --- |
| `0` | 4 | VRT header |
| `4` | 4 | Stream ID |
| `8` | 8 | Class ID |
| `16` | 4 | UTC integer seconds |
| `20` | 8 | fractional picoseconds |
| `28` | 4 | CIF0 |
| `32` | 4 | state indicators |
| `36` | 8 | payload format |
| `44` | 8 | sample rate |
| `52` | 8 | reference frequency |
| `60` | 8 | IF offset |
| `68` | 8 | bandwidth |
| `76` | 8 | ADC full-scale |
| `84` | 8 | receiver ID |
| `92` | variable | NUL-terminated, 32-bit-padded ASCII JSON metadata |

The serializer accepts exactly `kFersContextCif0` and rejects any other CIF0 value. Current CIF0 includes:

- State indicators
- Payload format
- Sample rate
- Reference frequency
- IF offset
- Bandwidth
- Reference level
- Device identifier
- ASCII metadata

Current payload format word:

```text
0x2110100200010001
```

The source comment defines that internal FERS word as:

```text
[63:60] complex Cartesian
[59:56] two's-complement integer
[55:48] item bits = 16
[47:40] packing bits = 16
[39:32] vector size = 2
[31:16] repeat count = 1
[15:0] reserved/version = 1
```

Context state indicators use the same indicator bit definitions as signal trailers. Context flags inside the ASCII metadata are FERS-specific JSON fields, not CIF0 bits.

Context flag bits:

| Bit | Name | Set when |
| ---: | --- | --- |
| `0` | `ContextFlagDechirped` | Receiver descriptor says dechirp is enabled. |
| `1` | `ContextFlagIfResampled` | Receiver descriptor says FMCW IF resampling is active. |
| `2` | `ContextFlagSampleLoss` | The stream has pending sample loss from a send failure. |
| `3` | `ContextFlagStreamOpen` | This context packet marks stream open. |
| `4` | `ContextFlagStreamClose` | This context packet marks stream close. |
| `5` | `ContextFlagFmcwMetadataPresent` | FMCW metadata is present and allowed by receiver mode. |
| `6` | `ContextFlagCwMetadataPresent` | CW metadata is present and allowed by receiver mode. |
| `7` | `ContextFlagPulsedMetadataPresent` | Pulsed metadata is present and allowed by receiver mode. |
| `8` | `ContextFlagSfcwMetadataPresent` | SFCW metadata is present and allowed by receiver mode. |

## Context ASCII JSON

The context ASCII metadata JSON always includes:

- `schema`: currently `fers-vita49-context-v1`
- `simulation_name`
- `receiver`
- `coordinate_frame`
- `initial_platform_state`
- `waveform`

`receiver` contains:

- `id`
- `name`
- `mode`
- `adc_bits`
- `context_flags`

`coordinate_frame` contains:

- `frame`
- `origin.latitude`
- `origin.longitude`
- `origin.altitude`
- `utm_zone`
- `utm_north_hemisphere`

`initial_platform_state` contains:

- `platform_id`
- `platform_name`
- `position_m.x`
- `position_m.y`
- `position_m.z`
- `velocity_mps.x`
- `velocity_mps.y`
- `velocity_mps.z`
- `rotation_rad.azimuth`
- `rotation_rad.elevation`

`waveform` contains a `kind` field. For known receiver modes it also contains `metadata_ref` pointing to the mode-specific object:

- `pulsed`
- `cw`
- `fmcw`
- `sfcw`

The serializer follows the receiver mode when stale mode blocks are present. For example, if receiver mode is `cw`, CW metadata is emitted and stale FMCW metadata is not emitted.

### Pulsed Metadata

`pulsed` contains:

- `present`
- `waveform_id`
- `waveform_name`
- `carrier_hz`
- `power_w`
- `pulse_width_s`
- `native_sample_rate_hz`
- `native_sample_count`
- `window_length_s`
- `window_prf_hz`
- `window_skip_s`
- `window_count`
- `pri_s`

`pri_s` is `1.0 / window_prf_hz` when PRF is positive; otherwise it is `null`.

### CW Metadata

`cw` contains:

- `present`
- `waveform_id`
- `waveform_name`
- `carrier_hz`
- `power_w`

### FMCW Metadata

`fmcw` contains:

- `present`
- `waveform_shape`
- `chirp_bandwidth_hz`
- `chirp_duration_s`
- `chirp_period_s`
- `chirp_rate_hz_per_s`
- `chirp_rate_signed_hz_per_s`
- `sweep_direction`
- `start_frequency_offset_hz`
- `triangle_period_s`
- `chirp_count`
- `triangle_count`
- `dechirp_mode`
- `dechirp_reference_source`
- `dechirp_reference_transmitter_id`
- `dechirp_reference_transmitter_name`
- `dechirp_reference_waveform_id`
- `dechirp_reference_waveform_name`

### SFCW Metadata

`sfcw` contains:

- `present`
- `waveform_id`
- `waveform_name`
- `carrier_hz`
- `start_frequency_offset_hz`
- `step_size_hz`
- `step_count`
- `dwell_time_s`
- `step_period_s`
- `sweep_period_s`
- `sweep_count`
- `first_frequency_hz`
- `last_frequency_hz`
- `frequency_span_hz`
- `effective_bandwidth_hz`

## Final Output Metadata

After a VITA run, `fers_get_last_output_metadata_json` includes a top-level `vita49` object. The object is omitted after HDF5 runs and after switching a context back to HDF5.

`vita49` fields:

- `endpoint`: `host:port`
- `endpoint_host`
- `endpoint_port`
- `epoch_unix_nanoseconds`: string when known, otherwise `null`
- `class_id`: currently `0xFA52530001000101`
- `adc_fullscale`
- `max_udp_payload`
- `queue_depth`
- `streams`

Each `streams` entry contains:

- `receiver_id`
- `receiver_name`
- `stream_id`
- `mode`
- `sample_rate`
- `reference_frequency`
- `packets_emitted`
- `samples_emitted`
- `packets_dropped`
- `samples_dropped`
- `over_range_count`
- `late_packet_count`
- `context_packet_count`
- `first_sample_time`
- `end_sample_time`
- `first_timestamp`
- `end_timestamp`

`first_timestamp` and `end_timestamp` are objects with:

- `integer_seconds`
- `fractional_picoseconds`

`end_sample_time` and `end_timestamp` are exclusive stream-end values from the current stats model.

## Live Telemetry Through C API

`fers_run_simulation_ex` can receive a VITA telemetry callback. The callback gets:

- Optional stream counter JSON.
- Optional packet trace batch JSON.

Stream counter JSON uses `context_packets` for live stats. Final output metadata uses `context_packet_count`.

Packet trace fields:

- `sequence`
- `event`: current source emits `data`, `context`, or `drop`
- `stream_id`
- `byte_count`
- `sample_count`
- `first_sample_time`
- `timestamp`
- `data_packet`
- `context_packet`
- `dropped`
- `over_range`
- `sample_loss`

The sink batches packet traces. Current constants are:

- Stats emit interval: `250 ms`
- Packet trace emit interval: `250 ms`
- Packet trace batch size: `64`

Disabling packet trace telemetry stops per-packet diagnostic record creation, but stream counter telemetry can still be emitted.

## Loss And Late Packet Counters

The sender records loss only when `DatagramSender::send(...)` throws or fails to send the full datagram. Queue fullness blocks the producer and does not count as loss.

When a send failure occurs:

- Data packet failures increment dropped data packet and dropped sample counters.
- Context packet failures increment dropped context packet counters internally.
- The sink marks sample loss pending for that stream.
- A later data or context packet carries the sample-loss indicator.

`late_packet_count` increments when a packet is sent more than 1 ms after its scheduled due time.
