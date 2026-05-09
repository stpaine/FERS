# XML Schema Reference

FERS scenarios are XML files with the `.fersxml` extension. A scenario describes the simulation time span, radar waveforms, antennas, timing sources, moving platforms, transmitters, receivers, monostatic radars, and targets.

This page documents the current user-facing XML format.

## Important Rules

- Element order matters. The order shown in this page is the order expected by validation.
- Most numeric values are written as text. Use normal decimal or scientific notation, such as `1000`, `1.0e6`, or `2.99792458e8`.
- Times are in seconds.
- Frequencies and sample rates are in Hz.
- Power is in watts.
- Distances and platform positions are always linear meter-valued coordinates. The signal simulation uses the numeric platform vectors directly; `<coordinatesystem>` is for KML/geospatial export, not for transforming simulation geometry.
- Angles are controlled by `<rotationangleunit>` for platform rotation values. Antenna-pattern XML axes declare their own angle unit.
- Names are used for references. Keep names unique and descriptive.
- XML schema validation checks structure. FERS also checks references, numeric limits, mode compatibility, and some physical constraints while loading or running.

## File Skeleton

```xml
<simulation name="ExampleScenario">
    <parameters>
        <starttime>0</starttime>
        <endtime>1.0</endtime>
        <rate>1000000</rate>
    </parameters>

    <waveform name="WaveformName">...</waveform>
    <timing name="TimingName">...</timing>
    <antenna name="AntennaName" pattern="isotropic"/>
    <platform name="PlatformName">...</platform>
</simulation>
```

`<parameters>` must come first. After that, waveforms, timings, antennas, platforms, and includes can appear in any order.

## `<simulation>`

Root element for one scenario.

| Attribute | Required | Meaning |
| --- | --- | --- |
| `name` | Yes | Human-readable scenario name. Used in logs, serialized output, and KML defaults. |

Children:

1. One `<parameters>` block.
2. Any number of `<waveform>`, `<timing>`, `<antenna>`, `<platform>`, and `<include>` elements.

A useful scenario normally needs at least one waveform, timing source, antenna, and platform.

## `<parameters>`

Global simulation settings.

Required children, in order:

| Element | Unit | Meaning |
| --- | --- | --- |
| `<starttime>` | seconds | Start time of the simulation. Usually `0`. |
| `<endtime>` | seconds | End time of the simulation. Must be after `starttime`. |
| `<rate>` | Hz | Main output sample rate. This is the sample rate used for final receiver data unless an FMCW IF chain sets `if_sample_rate`. |

Optional children, in order:

| Element | Unit | Default | Meaning |
| --- | --- | --- | --- |
| `<c>` | m/s | `299792458` | Propagation speed. Use the default for free-space radio propagation. |
| `<simSamplingRate>` | Hz | `1000` | Geometry interpolation rate used while building pulsed responses. Increase it when fast platform motion needs finer geometry sampling. |
| `<randomseed>` | integer | random | Seed for repeatable noise and random models. Set this when you want repeatable runs. |
| `<adc_bits>` | bits | `0` | ADC quantization depth. `0` disables quantization and stores normalized floating-point samples. |
| `<oversample>` | multiplier | `1` | Internal oversampling factor. Valid range is `1` to `8`. It does not normally change final output sample rate. |
| `<rotationangleunit>` | `deg` or `rad` | `deg` | Unit used by platform rotation values. |
| `<origin>` | attributes | built-in default | Geodetic reference used by ENU KML/geospatial export. |
| `<coordinatesystem>` | attributes | `ENU` | Coordinate frame used when converting platform coordinates to geodetic coordinates for KML/geospatial export. It does not transform positions during simulation. |

### `<origin>`

```xml
<origin latitude="-33.957652" longitude="18.4611991" altitude="111.01"/>
```

| Attribute | Required | Unit | Meaning |
| --- | --- | --- | --- |
| `latitude` | Yes | degrees | Geodetic latitude. |
| `longitude` | Yes | degrees | Geodetic longitude. |
| `altitude` | No | meters | Height above the ellipsoid or map reference. Defaults to `0` when an origin is supplied without altitude. |

If `<origin>` is omitted, FERS uses a built-in default origin. For real geospatial scenarios, set it explicitly.

### `<coordinatesystem>`

`<coordinatesystem>` controls how FERS converts platform coordinates to latitude, longitude, and altitude for KML/geospatial export. It does not convert waypoints into another internal frame before simulation. During simulation, FERS uses the numeric `x`, `y`, and `altitude` values directly for platform geometry and distance calculations.

```xml
<coordinatesystem frame="ENU"/>
```

or:

```xml
<coordinatesystem frame="UTM" zone="34" hemisphere="S"/>
```

| Attribute | Required | Meaning |
| --- | --- | --- |
| `frame` | Yes | Coordinate frame. Use `ENU`, `UTM`, or `ECEF`. |
| `zone` | Required for UTM | UTM zone, `1` to `60`. |
| `hemisphere` | Required for UTM | `N` or `S`. |

Frames:

| Frame | Meaning |
| --- | --- |
| `ENU` | KML export treats platform `x`, `y`, and `altitude` as east, north, and up offsets in meters from `<origin>`. This is the easiest frame for local scenes. |
| `UTM` | KML export treats platform `x` and `y` as UTM easting and northing in meters. `altitude` is also meters. Set `zone` and `hemisphere`. |
| `ECEF` | KML export treats platform `x`, `y`, and `altitude` as Earth-centered, Earth-fixed Cartesian coordinates in meters. In this frame, the XML element name `altitude` represents ECEF Z, not geodetic altitude. |

FERS does not accept platform latitude/longitude waypoints. Latitude and longitude appear only in `<origin>`, where they define the geodetic reference used by ENU and KML export.

## `<waveform>`

A waveform defines what a transmitter emits.

```xml
<waveform name="PulseWave">
    <power>10000</power>
    <carrier_frequency>10.0e9</carrier_frequency>
    <pulsed_from_file filename="pulse.h5"/>
</waveform>
```

| Attribute | Required | Meaning |
| --- | --- | --- |
| `name` | Yes | Name used by transmitters and monostatic radars. |

Required children:

| Element | Unit | Meaning |
| --- | --- | --- |
| `<power>` | watts | Transmit power. |
| `<carrier_frequency>` | Hz | RF carrier frequency. The waveform data itself is complex baseband. |

Then choose exactly one waveform type:

- `<pulsed_from_file>`.
- `<cw/>`.
- `<fmcw_linear_chirp>`.
- `<fmcw_triangle>`.

The chosen waveform type must match the radar mode that uses it.

### `<pulsed_from_file>`

Use this for pulsed radar.

```xml
<pulsed_from_file filename="pulse.h5"/>
```

| Attribute | Required | Meaning |
| --- | --- | --- |
| `filename` | Yes | Pulse waveform file. Relative paths are resolved relative to the scenario file. |

Supported file types:

| Type | Expected content |
| --- | --- |
| HDF5 `.h5` | Datasets `/I/value` and `/Q/value`. |
| CSV `.csv` | Sample count, sample rate, then complex samples. |

The file contains complex baseband samples. Do not include the RF carrier in the file; use `<carrier_frequency>` for that.

### `<cw/>`

Use this for continuous-wave simulation.

```xml
<cw/>
```

There are no child parameters. The signal is a continuous tone at the chosen carrier frequency.

### `<fmcw_linear_chirp>`

Use this for sawtooth up-chirp or down-chirp FMCW.

```xml
<fmcw_linear_chirp direction="up">
    <chirp_bandwidth>20.0e6</chirp_bandwidth>
    <chirp_duration>1.0e-3</chirp_duration>
    <chirp_period>1.0e-3</chirp_period>
    <start_frequency_offset>0</start_frequency_offset>
    <chirp_count>64</chirp_count>
</fmcw_linear_chirp>
```

| Attribute | Required | Meaning |
| --- | --- | --- |
| `direction` | Yes | `up` increases frequency through the chirp. `down` decreases frequency. |

| Element | Required | Unit | Meaning |
| --- | --- | --- | --- |
| `<chirp_bandwidth>` | Yes | Hz | Frequency sweep width. Must be positive. |
| `<chirp_duration>` | Yes | seconds | Active sweep duration. Must be positive. |
| `<chirp_period>` | Yes | seconds | Time from one chirp start to the next. Must be at least `chirp_duration`. |
| `<start_frequency_offset>` | No | Hz | Baseband offset at the start of the chirp. Defaults to `0`. |
| `<chirp_count>` | No | count | Maximum chirps per schedule period. If omitted, the schedule duration determines how many chirps fit. |

FMCW chirps must fit the selected sample rate. If the sweep bandwidth is too large for `<rate> * <oversample>`, FERS rejects or warns about the scenario.

### `<fmcw_triangle>`

Use this for triangular FMCW sweeps.

```xml
<fmcw_triangle>
    <chirp_bandwidth>20.0e6</chirp_bandwidth>
    <chirp_duration>1.0e-3</chirp_duration>
    <triangle_count>64</triangle_count>
</fmcw_triangle>
```

| Element | Required | Unit | Meaning |
| --- | --- | --- | --- |
| `<chirp_bandwidth>` | Yes | Hz | Sweep width for each leg. |
| `<chirp_duration>` | Yes | seconds | Duration of one leg. A full triangle takes `2 * chirp_duration`. |
| `<start_frequency_offset>` | No | Hz | Baseband offset at the start of the triangle. Defaults to `0`. |
| `<triangle_count>` | No | count | Maximum complete triangles per schedule period. |

## `<timing>`

A timing source represents the clock or oscillator used by transmitters and receivers.

```xml
<timing name="Clock" synconpulse="false">
    <frequency>10.0e6</frequency>
    <freq_offset>0</freq_offset>
    <phase_offset>0</phase_offset>
</timing>
```

| Attribute | Required | Default | Meaning |
| --- | --- | --- | --- |
| `name` | Yes | none | Name used by transmitters, receivers, and monostatic radars. |
| `synconpulse` | No | `false` | For pulsed operation, controls whether timing is synchronized on pulse events. Use `false` unless you specifically need pulse-synchronized timing behavior. |

Children, in order:

| Element | Required | Unit | Meaning |
| --- | --- | --- | --- |
| `<frequency>` | Yes | Hz | Nominal oscillator frequency. |
| `<freq_offset>` | No | Hz | Constant frequency offset. |
| `<random_freq_offset_stdev>` | No | Hz | Standard deviation for random frequency offset. |
| `<phase_offset>` | No | radians | Constant phase offset. |
| `<random_phase_offset_stdev>` | No | radians | Standard deviation for random phase offset. |
| `<noise_entry>` | No | see below | Phase-noise model entry. Multiple entries are allowed. |

### `<noise_entry>`

```xml
<noise_entry>
    <alpha>0</alpha>
    <weight>-100</weight>
</noise_entry>
```

| Element | Meaning |
| --- | --- |
| `<alpha>` | Phase-noise slope category used by the timing model. |
| `<weight>` | Noise weighting for that category, conventionally written in dBc/Hz-style values. |

Use timing noise only when your analysis needs oscillator effects. For first scenarios, leave the optional timing fields out.

## `<antenna>`

An antenna defines directional gain.

```xml
<antenna name="Iso" pattern="isotropic">
    <efficiency>1.0</efficiency>
</antenna>
```

| Attribute | Required | Meaning |
| --- | --- | --- |
| `name` | Yes | Name used by transmitters, receivers, and monostatic radars. |
| `pattern` | Yes | Pattern model. Supported values are `isotropic`, `sinc`, `gaussian`, `squarehorn`, `parabolic`, `xml`, and `file`. |
| `filename` | Depends | Required by `xml` and `file` patterns. |

Optional child parameters:

| Element | Used by | Meaning |
| --- | --- | --- |
| `<alpha>` | `sinc` | Gain scale for the sinc pattern. |
| `<beta>` | `sinc` | Angular scale for the sinc pattern. |
| `<gamma>` | `sinc` | Shape exponent for the sinc pattern. |
| `<diameter>` | `squarehorn`, `parabolic` | Aperture diameter or size in meters. |
| `<azscale>` | `gaussian` | Azimuth spread for Gaussian pattern. |
| `<elscale>` | `gaussian` | Elevation spread for Gaussian pattern. |
| `<efficiency>` | most patterns | Efficiency multiplier. Defaults to `1`. |

Pattern notes:

| Pattern | Meaning |
| --- | --- |
| `isotropic` | Equal gain in every direction. Good for simple scenarios and debugging. |
| `sinc` | Analytic sinc-shaped gain pattern controlled by `alpha`, `beta`, and `gamma`. |
| `gaussian` | Analytic Gaussian beam controlled by `azscale` and `elscale`. |
| `squarehorn` | Approximate square horn aperture model using `diameter`. |
| `parabolic` | Approximate parabolic dish model using `diameter`. |
| `xml` | Uses a standalone antenna-pattern XML file. |
| `file` | Uses an HDF5 antenna pattern file. |

## Standalone Antenna Pattern XML

Use `pattern="xml"` when you want to provide azimuth and elevation gain samples.

```xml
<antenna>
    <azimuth unit="deg" format="dBi" symmetry="mirrored">
        <gainsample>
            <angle>0</angle>
            <gain>20</gain>
        </gainsample>
    </azimuth>
    <elevation unit="deg" format="dBi" symmetry="mirrored">
        <gainsample>
            <angle>0</angle>
            <gain>20</gain>
        </gainsample>
    </elevation>
</antenna>
```

Axis attributes:

| Attribute | Values | Default | Meaning |
| --- | --- | --- | --- |
| `unit` | `rad`, `deg` | `rad` | Angle unit for this axis file. |
| `format` | `linear`, `dBi` | `linear` | Gain value format. |
| `symmetry` | `mirrored`, `none` | `mirrored` | Whether positive-angle samples are mirrored to negative angles. |

Each `<gainsample>` contains:

| Element | Meaning |
| --- | --- |
| `<angle>` | Angle on that axis. |
| `<gain>` | Gain at that angle. |

## `<platform>`

A platform is anything that has position and orientation: radar vehicle, receiver site, aircraft, target, or any other moving object.

```xml
<platform name="RadarPlatform">
    <motionpath interpolation="static">
        <positionwaypoint>
            <x>0</x>
            <y>0</y>
            <altitude>0</altitude>
            <time>0</time>
        </positionwaypoint>
    </motionpath>
    <fixedrotation>
        <startazimuth>0</startazimuth>
        <startelevation>0</startelevation>
        <azimuthrate>0</azimuthrate>
        <elevationrate>0</elevationrate>
    </fixedrotation>
</platform>
```

| Attribute | Required | Meaning |
| --- | --- | --- |
| `name` | Yes | Platform name. |

Children, in order:

1. One `<motionpath>`.
2. Either `<rotationpath>` or `<fixedrotation>`.
3. Zero or more `<monostatic>`, `<transmitter>`, `<receiver>`, or `<target>` elements.

### `<motionpath>`

| Attribute | Values | Default | Meaning |
| --- | --- | --- | --- |
| `interpolation` | `static`, `linear`, `cubic` | `static` | How FERS interpolates between position waypoints. |

Contains one or more `<positionwaypoint>` elements.

### `<positionwaypoint>`

| Element | Unit | Meaning |
| --- | --- | --- |
| `<x>` | meters | X coordinate used directly by the simulation. KML export treats it as ENU east offset, UTM easting, or ECEF X depending on `<coordinatesystem>`. |
| `<y>` | meters | Y coordinate used directly by the simulation. KML export treats it as ENU north offset, UTM northing, or ECEF Y depending on `<coordinatesystem>`. |
| `<altitude>` | meters | Z coordinate used directly by the simulation. KML export treats it as ENU up offset, UTM altitude, or ECEF Z depending on `<coordinatesystem>`. |
| `<time>` | seconds | Time when the platform is at this waypoint. |

For `static`, the first waypoint is normally enough. For `linear` or `cubic`, provide enough waypoints to describe the path over the simulation time.

### `<rotationpath>`

Use rotation waypoints when pointing changes according to specific time points.

| Attribute | Values | Meaning |
| --- | --- | --- |
| `interpolation` | `static`, `linear`, `cubic` | How FERS interpolates between rotation waypoints. |

Contains one or more `<rotationwaypoint>` elements.

| Element | Unit | Meaning |
| --- | --- | --- |
| `<azimuth>` | `<rotationangleunit>` | Platform azimuth. |
| `<elevation>` | `<rotationangleunit>` | Platform elevation. |
| `<time>` | seconds | Time for this pointing waypoint. |

### `<fixedrotation>`

Use fixed rotation when pointing changes at a constant rate.

| Element | Unit | Meaning |
| --- | --- | --- |
| `<startazimuth>` | `<rotationangleunit>` | Azimuth at simulation time zero. |
| `<startelevation>` | `<rotationangleunit>` | Elevation at simulation time zero. |
| `<azimuthrate>` | `<rotationangleunit>` per second | Constant azimuth rate. |
| `<elevationrate>` | `<rotationangleunit>` per second | Constant elevation rate. |

## `<schedule>`

Schedules are optional. They restrict when a transmitter or receiver is active.

```xml
<schedule>
    <period start="0.0" end="0.5"/>
    <period start="0.8" end="1.0"/>
</schedule>
```

| Element or attribute | Unit | Meaning |
| --- | --- | --- |
| `<period>` | none | One active time interval. |
| `start` | seconds | Start time of the active interval. |
| `end` | seconds | End time of the active interval. Must be after `start`. |

If no schedule is supplied, the object is active for the full simulation.

## Radar Modes

Each transmitter, receiver, or monostatic radar chooses exactly one mode:

- `<pulsed_mode>`.
- `<cw_mode/>`.
- `<fmcw_mode>`.

The mode must match the waveform type.

## `<transmitter>`

Use a transmitter for bistatic or multistatic layouts.

```xml
<transmitter name="Tx" waveform="Wave" antenna="TxAntenna" timing="TxClock">
    <pulsed_mode>
        <prf>1000</prf>
    </pulsed_mode>
</transmitter>
```

| Attribute | Required | Meaning |
| --- | --- | --- |
| `name` | Yes | Transmitter name. |
| `waveform` | Yes | Name of a `<waveform>`. |
| `antenna` | Yes | Name of an `<antenna>`. |
| `timing` | Yes | Name of a `<timing>` source. |

Children:

1. Exactly one mode block.
2. Optional `<schedule>`.

### Transmitter `<pulsed_mode>`

| Element | Unit | Meaning |
| --- | --- | --- |
| `<prf>` | Hz | Pulse repetition frequency. Use a positive value. |

### Transmitter `<cw_mode/>`

No child parameters.

### Transmitter `<fmcw_mode>`

No child parameters for transmitters. Dechirp and IF settings belong on receivers or monostatic radars.

## `<receiver>`

Use a receiver for bistatic or multistatic layouts.

```xml
<receiver name="Rx" antenna="RxAntenna" timing="RxClock" nodirect="false" nopropagationloss="false">
    <pulsed_mode>
        <prf>1000</prf>
        <window_skip>0</window_skip>
        <window_length>0.001</window_length>
    </pulsed_mode>
    <noise_temp>290</noise_temp>
</receiver>
```

| Attribute | Required | Default | Meaning |
| --- | --- | --- | --- |
| `name` | Yes | none | Receiver name. Also used in output file names. |
| `antenna` | Yes | none | Name of an `<antenna>`. |
| `timing` | Yes | none | Name of a `<timing>` source. |
| `nodirect` | No | `false` | If `true`, direct transmitter-to-receiver paths are ignored. Target reflections can still be received. |
| `nopropagationloss` | No | `false` | If `true`, free-space propagation-loss scaling is disabled for this receiver. Useful for debugging, not realistic link budgets. |

Children:

1. Exactly one mode block.
2. Optional `<noise_temp>`.
3. Optional `<schedule>`.

### Receiver `<pulsed_mode>`

| Element | Unit | Meaning |
| --- | --- | --- |
| `<prf>` | Hz | Receiver window repetition frequency. Must be positive. |
| `<window_skip>` | seconds | Delay from each window trigger before recording starts. Must be zero or positive. |
| `<window_length>` | seconds | Recording duration for each receiver window. Must be positive. |

Use `window_skip` to choose the range region of interest.

### Receiver `<cw_mode/>`

No child parameters.

### Receiver `<fmcw_mode>`

See the "FMCW Mode: Receiver And Monostatic" section below.

### `<noise_temp>`

| Element | Unit | Default | Meaning |
| --- | --- | --- | --- |
| `<noise_temp>` | kelvin | `0` | Receiver noise temperature. Use `290` K as a common room-temperature reference. |

## `<monostatic>`

Use a monostatic radar when transmitter and receiver are colocated on the same platform.

```xml
<monostatic name="Radar" antenna="Antenna" waveform="Wave" timing="Clock">
    <pulsed_mode>
        <prf>1000</prf>
        <window_skip>0</window_skip>
        <window_length>0.001</window_length>
    </pulsed_mode>
    <noise_temp>290</noise_temp>
</monostatic>
```

| Attribute | Required | Default | Meaning |
| --- | --- | --- | --- |
| `name` | Yes | none | Radar name. Also used in output file names. |
| `antenna` | Yes | none | Shared antenna. |
| `waveform` | Yes | none | Transmitted waveform. |
| `timing` | Yes | none | Shared timing source. |
| `nodirect` | No | `false` | If `true`, the direct self path is ignored. |
| `nopropagationloss` | No | `false` | If `true`, propagation-loss scaling is disabled. Useful for debugging, not realistic link budgets. |

Children:

1. Exactly one mode block.
2. Optional `<noise_temp>`.
3. Optional `<schedule>`.

Mode blocks are the same as receiver mode blocks, except that a monostatic pulsed mode also supplies the transmitter PRF.

## FMCW Mode: Receiver And Monostatic

FMCW receivers can either export the raw received signal or dechirp it inside FERS.

```xml
<fmcw_mode dechirp_mode="physical">
    <dechirp_reference source="attached"/>
    <if_sample_rate>2000000</if_sample_rate>
    <if_filter_bandwidth>800000</if_filter_bandwidth>
    <if_filter_transition_width>100000</if_filter_transition_width>
</fmcw_mode>
```

| Attribute | Values | Default | Meaning |
| --- | --- | --- | --- |
| `dechirp_mode` | `none`, `physical`, `ideal` | `none` | Controls whether and how received FMCW is mixed down to IF. |

Modes:

| Mode | Meaning |
| --- | --- |
| `none` | Do not dechirp. Output raw received complex baseband. Do not provide `dechirp_reference` or IF settings. |
| `physical` | Dechirp using a reference signal while preserving timing effects between source and receiver. Use for more realistic IF output. |
| `ideal` | Dechirp with an ideal reference. Useful when you want the beat signal without timing-source decorrelation. |

### `<dechirp_reference>`

Required for `physical` and `ideal`. Not allowed for `none`.

| Attribute | Required | Meaning |
| --- | --- | --- |
| `source` | Yes | `attached`, `transmitter`, or `custom`. |
| `transmitter_name` | For `source="transmitter"` | Name of the transmitter whose waveform is used as the reference. |
| `waveform_name` | For `source="custom"` | Name of the waveform used as the reference. |

Reference sources:

| Source | Meaning |
| --- | --- |
| `attached` | Use the waveform attached to the monostatic radar. This is the normal monostatic choice. |
| `transmitter` | Use a named transmitter as the reference. Useful for bistatic scenes. |
| `custom` | Use a named waveform directly as the reference. |

### IF Settings

| Element | Unit | Meaning |
| --- | --- | --- |
| `<if_sample_rate>` | Hz | Optional output sample rate after dechirp. Must be positive and no higher than the main simulation sample rate. |
| `<if_filter_bandwidth>` | Hz | Optional IF low-pass bandwidth. Requires `if_sample_rate` and must be less than half of it. |
| `<if_filter_transition_width>` | Hz | Optional transition width for the IF filter. Requires `if_sample_rate`. |

If you omit `if_filter_bandwidth` but request IF resampling, FERS chooses a conservative default based on the output rate.

## `<target>`

A target is a radar scatterer attached to a platform.

```xml
<target name="Target">
    <rcs type="isotropic">
        <value>10</value>
    </rcs>
</target>
```

| Attribute | Required | Meaning |
| --- | --- | --- |
| `name` | Yes | Target name. |

Children:

1. Required `<rcs>`.
2. Optional `<model>`.

### `<rcs>`

| Attribute | Required | Meaning |
| --- | --- | --- |
| `type` | Yes | `isotropic` or `file`. |
| `filename` | For `type="file"` | Target RCS pattern file. |

For constant RCS:

```xml
<rcs type="isotropic">
    <value>10</value>
</rcs>
```

| Element | Unit | Meaning |
| --- | --- | --- |
| `<value>` | square meters | Constant radar cross section. |

For file-backed RCS:

```xml
<rcs type="file" filename="target_rcs.xml"/>
```

Target RCS files use azimuth and elevation sample axes with `<rcssample>` entries containing `<angle>` and `<rcs>`.

### `<model>`

Optional target fluctuation model.

```xml
<model type="chisquare">
    <k>2</k>
</model>
```

| Attribute | Values | Meaning |
| --- | --- | --- |
| `type` | `constant`, `chisquare`, `gamma` | RCS fluctuation model. |

| Element | Used by | Meaning |
| --- | --- | --- |
| `<k>` | `chisquare`, `gamma` | Shape or degrees-of-freedom parameter for the fluctuation model. |

If `<model>` is omitted, RCS is constant.

## `<include>`

Use includes to split large scenarios into separate files.

```xml
<include>assets.xml</include>
<include>platforms.xml</include>
```

The included file path is resolved relative to the main scenario file. Includes are expanded before validation when loading from a file.

Includes are not expanded when loading an XML string through `libfers`.

## Practical Constraints Not Fully Expressed By XML Schema

The XML schema is intentionally simple, so some important checks happen when FERS loads or runs the scenario:

- Names referenced by `waveform`, `antenna`, and `timing` must exist.
- A waveform can be used only with a compatible mode.
- `<rate>` must be positive.
- `<endtime>` must be after `<starttime>`.
- `<oversample>` must be between `1` and `8`.
- Pulsed receiver `prf` and `window_length` must be positive.
- Pulsed receiver `window_skip` must not be negative.
- FMCW chirp bandwidth and duration must be positive.
- FMCW chirp period must be at least chirp duration.
- FMCW sweep bandwidth must fit the effective sample rate.
- Receiver noise temperature must not be negative.
- UTM coordinates need a valid zone and hemisphere.
- File-backed waveform, antenna, and RCS assets must exist and match the expected file structure.
