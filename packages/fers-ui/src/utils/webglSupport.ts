// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

export type WebGLContextMode = 'webgl2' | 'webgl' | 'experimental-webgl';

interface NavigatorWithUserAgentData extends Navigator {
    userAgentData?: {
        platform?: string;
    };
}

export interface WebGLProbeResult {
    mode: WebGLContextMode;
    available: boolean;
    usable: boolean;
    lostAtCreation: boolean;
    lostAfterExercise: boolean;
    clearError: number | null;
    version: string | null;
    renderer: string | null;
    vendor: string | null;
    shaderCreated: boolean;
    reason: string;
}

export interface WebGLPlatformInfo {
    platform: string;
    userAgent: string;
    isMac: boolean;
    isIntelMac: boolean;
}

export interface WebGLSupportReport {
    rendererSupported: boolean;
    supportedMode: WebGLContextMode | null;
    summary: string;
    platform: WebGLPlatformInfo;
    probes: WebGLProbeResult[];
}

type WebGLProbeContext = WebGLRenderingContext | WebGL2RenderingContext;

const CONTEXT_MODES: readonly WebGLContextMode[] = [
    'webgl2',
    'webgl',
    'experimental-webgl',
];

let cachedReportPromise: Promise<WebGLSupportReport> | null = null;

function getPlatformInfo(): WebGLPlatformInfo {
    if (typeof navigator === 'undefined') {
        return {
            platform: 'unknown',
            userAgent: 'unknown',
            isMac: false,
            isIntelMac: false,
        };
    }

    const nav = navigator as NavigatorWithUserAgentData;
    const platform =
        nav.userAgentData?.platform ?? navigator.platform ?? 'unknown';
    const userAgent = navigator.userAgent ?? 'unknown';
    const platformText = `${platform} ${userAgent}`;
    const isMac = /mac/i.test(platformText);
    const isIntelMac = isMac && /(intel|x86_64|macintel)/i.test(platformText);

    return {
        platform,
        userAgent,
        isMac,
        isIntelMac,
    };
}

function safeIsContextLost(gl: WebGLProbeContext): boolean {
    try {
        return typeof gl.isContextLost === 'function'
            ? gl.isContextLost()
            : false;
    } catch {
        return true;
    }
}

function safeGetParameter(
    gl: WebGLProbeContext,
    parameter: number
): string | null {
    try {
        const value = gl.getParameter(parameter);

        if (value === null || value === undefined) {
            return null;
        }

        return String(value);
    } catch {
        return null;
    }
}

function releaseContext(gl: WebGLProbeContext): void {
    try {
        gl.getExtension('WEBGL_lose_context')?.loseContext();
    } catch {
        // Best-effort cleanup only.
    }
}

function buildFailureReason(result: {
    available: boolean;
    lostAtCreation: boolean;
    lostAfterExercise: boolean;
    clearError: number | null;
    shaderCreated: boolean;
    version: string | null;
}): string {
    if (!result.available) {
        return 'Context acquisition returned null.';
    }

    if (result.lostAtCreation) {
        return 'Context was already lost at creation time.';
    }

    if (result.clearError === 37442) {
        return 'Basic GL commands returned CONTEXT_LOST_WEBGL.';
    }

    if (result.lostAfterExercise) {
        return 'Context was lost during the startup probe.';
    }

    if (!result.shaderCreated) {
        return 'Shader creation failed during the startup probe.';
    }

    if (result.version === null) {
        return 'Capability queries returned null.';
    }

    return 'Renderer startup checks failed.';
}

function probeContext(mode: WebGLContextMode): WebGLProbeResult {
    if (typeof document === 'undefined') {
        return {
            mode,
            available: false,
            usable: false,
            lostAtCreation: false,
            lostAfterExercise: false,
            clearError: null,
            version: null,
            renderer: null,
            vendor: null,
            shaderCreated: false,
            reason: 'Document is unavailable in the current runtime.',
        };
    }

    const canvas = document.createElement('canvas');
    canvas.width = 1;
    canvas.height = 1;

    let gl: WebGLProbeContext | null = null;

    try {
        gl = canvas.getContext(mode, {
            alpha: true,
            antialias: false,
            depth: true,
            stencil: false,
            powerPreference: 'default',
        }) as WebGLProbeContext | null;

        if (!gl) {
            return {
                mode,
                available: false,
                usable: false,
                lostAtCreation: false,
                lostAfterExercise: false,
                clearError: null,
                version: null,
                renderer: null,
                vendor: null,
                shaderCreated: false,
                reason: 'Context acquisition returned null.',
            };
        }

        const lostAtCreation = safeIsContextLost(gl);
        const version = safeGetParameter(gl, gl.VERSION);
        const renderer = safeGetParameter(gl, gl.RENDERER);
        const vendor = safeGetParameter(gl, gl.VENDOR);

        let clearError: number | null = null;
        let lostAfterExercise = lostAtCreation;
        let shaderCreated = false;

        try {
            gl.viewport(0, 0, 1, 1);
            gl.clearColor(0, 0, 0, 1);
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
            clearError = gl.getError();
            lostAfterExercise = safeIsContextLost(gl);

            const shader = gl.createShader(gl.VERTEX_SHADER);
            shaderCreated = shader !== null;

            if (shader) {
                gl.deleteShader(shader);
            }
        } catch {
            lostAfterExercise = safeIsContextLost(gl);
        }

        const usable =
            !lostAtCreation &&
            !lostAfterExercise &&
            clearError !== gl.CONTEXT_LOST_WEBGL &&
            shaderCreated &&
            version !== null;

        return {
            mode,
            available: true,
            usable,
            lostAtCreation,
            lostAfterExercise,
            clearError,
            version,
            renderer,
            vendor,
            shaderCreated,
            reason: usable
                ? 'Context passed startup checks.'
                : buildFailureReason({
                      available: true,
                      lostAtCreation,
                      lostAfterExercise,
                      clearError,
                      shaderCreated,
                      version,
                  }),
        };
    } catch (error) {
        return {
            mode,
            available: gl !== null,
            usable: false,
            lostAtCreation: gl ? safeIsContextLost(gl) : false,
            lostAfterExercise: gl ? safeIsContextLost(gl) : false,
            clearError: null,
            version: null,
            renderer: null,
            vendor: null,
            shaderCreated: false,
            reason:
                error instanceof Error
                    ? error.message
                    : 'Unknown WebGL startup error.',
        };
    } finally {
        if (gl) {
            releaseContext(gl);
        }

        canvas.width = 0;
        canvas.height = 0;
    }
}

function buildSummary(
    webgl2Probe: WebGLProbeResult | undefined,
    fallbackProbe: WebGLProbeResult | undefined,
    platform: WebGLPlatformInfo
): string {
    if (!webgl2Probe) {
        return 'WebGL2 startup probe did not run.';
    }

    if (!webgl2Probe.available) {
        return fallbackProbe?.usable
            ? 'Only legacy WebGL is usable on this system. The current 3D renderer requires a working WebGL2 context.'
            : 'No usable WebGL context was detected during startup.';
    }

    if (webgl2Probe.lostAtCreation || webgl2Probe.clearError === 37442) {
        return platform.isIntelMac
            ? 'The system WebKit/WebGL stack reported a lost WebGL context during startup on this Intel Mac.'
            : 'The system WebGL stack reported a lost WebGL context during startup.';
    }

    return webgl2Probe.reason;
}

function evaluateWebGLSupport(): WebGLSupportReport {
    const platform = getPlatformInfo();
    const probes = CONTEXT_MODES.map(probeContext);
    const webgl2Probe = probes.find((probe) => probe.mode === 'webgl2');
    const fallbackProbe = probes.find(
        (probe) => probe.mode !== 'webgl2' && probe.usable
    );

    if (webgl2Probe?.usable) {
        return {
            rendererSupported: true,
            supportedMode: 'webgl2',
            summary: 'WebGL2 startup probe succeeded.',
            platform,
            probes,
        };
    }

    return {
        rendererSupported: false,
        supportedMode: null,
        summary: buildSummary(webgl2Probe, fallbackProbe, platform),
        platform,
        probes,
    };
}

export function getWebGLSupportReport(): Promise<WebGLSupportReport> {
    if (!cachedReportPromise) {
        cachedReportPromise = Promise.resolve().then(evaluateWebGLSupport);
    }

    return cachedReportPromise;
}

export function resetWebGLSupportReportCache(): void {
    cachedReportPromise = null;
}
