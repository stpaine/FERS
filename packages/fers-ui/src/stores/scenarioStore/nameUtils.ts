// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import type { PlatformComponent, ScenarioData } from './types';

type ScenarioNameSource = Pick<
    ScenarioData,
    'waveforms' | 'timings' | 'antennas' | 'platforms'
>;

type UniqueNameOptions = {
    copy?: boolean;
};

export function getComponentIdentityIds(
    component: PlatformComponent
): string[] {
    return component.type === 'monostatic'
        ? [component.id, component.txId, component.rxId]
        : [component.id];
}

function toIgnoredSet(
    ignoredIds: Iterable<string | null | undefined> = []
): Set<string> {
    const ignored = new Set<string>();
    for (const id of ignoredIds) {
        if (id) {
            ignored.add(id);
        }
    }
    return ignored;
}

function shouldIgnore(
    ids: Iterable<string | null | undefined>,
    ignoredIds: Set<string>
): boolean {
    for (const id of ids) {
        if (id && ignoredIds.has(id)) {
            return true;
        }
    }
    return false;
}

export function collectScenarioNames(
    scenarioData: ScenarioNameSource,
    ignoredIds: Iterable<string | null | undefined> = []
): Set<string> {
    const ignored = toIgnoredSet(ignoredIds);
    const names = new Set<string>();

    const addName = (
        name: string,
        ids: Iterable<string | null | undefined>
    ) => {
        if (!shouldIgnore(ids, ignored)) {
            names.add(name);
        }
    };

    for (const waveform of scenarioData.waveforms) {
        addName(waveform.name, [waveform.id]);
    }
    for (const timing of scenarioData.timings) {
        addName(timing.name, [timing.id]);
    }
    for (const antenna of scenarioData.antennas) {
        addName(antenna.name, [antenna.id]);
    }
    for (const platform of scenarioData.platforms) {
        addName(platform.name, [platform.id]);
        for (const component of platform.components) {
            addName(component.name, getComponentIdentityIds(component));
        }
    }

    return names;
}

export function createUniqueName(
    baseName: string,
    existingNames: Set<string>,
    options: UniqueNameOptions = {}
): string {
    if (!options.copy && !existingNames.has(baseName)) {
        return baseName;
    }

    let candidate = `${baseName} Copy`;
    let suffix = 2;

    while (existingNames.has(candidate)) {
        candidate = `${baseName} Copy ${suffix}`;
        suffix += 1;
    }

    return candidate;
}

export function createUniqueScenarioName(
    scenarioData: ScenarioNameSource,
    baseName: string,
    ignoredIds: Iterable<string | null | undefined> = []
): string {
    return createUniqueName(
        baseName,
        collectScenarioNames(scenarioData, ignoredIds)
    );
}

export function createUniqueScenarioCopyName(
    scenarioData: ScenarioNameSource,
    baseName: string,
    ignoredIds: Iterable<string | null | undefined> = []
): string {
    return createUniqueName(
        baseName,
        collectScenarioNames(scenarioData, ignoredIds),
        { copy: true }
    );
}
