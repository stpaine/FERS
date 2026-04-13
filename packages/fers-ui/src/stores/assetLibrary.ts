// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

import { create } from 'zustand';
import {
    type AssetLibraryTemplate,
    createAssetLibraryFile,
    createTemplateFromScenarioItem,
    parseAssetTemplates,
    prepareTemplatesForCatalog,
} from './assetTemplates';
import { useScenarioStore } from './scenarioStore';

const CATALOG_FILENAME = 'asset-library.json';
const LOCAL_STORAGE_KEY = 'fers.assetLibrary.v1';

type AssetLibraryState = {
    templates: AssetLibraryTemplate[];
    isLoading: boolean;
    hasLoaded: boolean;
    error: string | null;
    loadCatalog: () => Promise<void>;
    saveScenarioItem: (
        itemId: string,
        nameOverride?: string
    ) => Promise<AssetLibraryTemplate | null>;
    addTemplates: (templates: AssetLibraryTemplate[]) => Promise<number>;
    updateTemplateName: (templateId: string, name: string) => Promise<void>;
    deleteTemplate: (templateId: string) => Promise<void>;
};

function canUseLocalStorage(): boolean {
    return typeof globalThis.localStorage !== 'undefined';
}

function readLocalCatalog(): AssetLibraryTemplate[] {
    if (!canUseLocalStorage()) {
        return [];
    }

    const raw = globalThis.localStorage.getItem(LOCAL_STORAGE_KEY);
    if (!raw) {
        return [];
    }

    return parseAssetTemplates(JSON.parse(raw));
}

function writeLocalCatalog(templates: AssetLibraryTemplate[]): void {
    if (!canUseLocalStorage()) {
        return;
    }

    globalThis.localStorage.setItem(
        LOCAL_STORAGE_KEY,
        JSON.stringify(createAssetLibraryFile(templates), null, 2)
    );
}

function isMissingCatalogError(error: unknown): boolean {
    const message = error instanceof Error ? error.message : String(error);
    return (
        message.includes('not found') ||
        message.includes('No such file') ||
        message.includes('os error 2') ||
        message.includes('ENOENT')
    );
}

async function getTauriCatalogLocation(): Promise<{
    directory: string;
    filePath: string;
}> {
    const { appDataDir, join } = await import('@tauri-apps/api/path');
    const directory = await appDataDir();
    return {
        directory,
        filePath: await join(directory, CATALOG_FILENAME),
    };
}

async function readCatalog(): Promise<AssetLibraryTemplate[]> {
    try {
        const { readTextFile } = await import('@tauri-apps/plugin-fs');
        const { filePath } = await getTauriCatalogLocation();
        const raw = await readTextFile(filePath);
        return parseAssetTemplates(JSON.parse(raw));
    } catch (error) {
        if (isMissingCatalogError(error)) {
            return [];
        }
        return readLocalCatalog();
    }
}

async function writeCatalog(templates: AssetLibraryTemplate[]): Promise<void> {
    try {
        const { mkdir, writeTextFile } = await import('@tauri-apps/plugin-fs');
        const { directory, filePath } = await getTauriCatalogLocation();
        await mkdir(directory, { recursive: true });
        await writeTextFile(
            filePath,
            JSON.stringify(createAssetLibraryFile(templates), null, 2)
        );
    } catch (_error) {
        writeLocalCatalog(templates);
    }
}

function sortTemplates(
    templates: AssetLibraryTemplate[]
): AssetLibraryTemplate[] {
    return [...templates].sort((a, b) => {
        const kindCompare = a.kind.localeCompare(b.kind);
        if (kindCompare !== 0) {
            return kindCompare;
        }
        return a.name.localeCompare(b.name);
    });
}

export const useAssetLibraryStore = create<AssetLibraryState>()((set, get) => ({
    templates: [],
    isLoading: false,
    hasLoaded: false,
    error: null,
    loadCatalog: async () => {
        set({ isLoading: true, error: null });
        try {
            const templates = await readCatalog();
            set({
                templates: sortTemplates(templates),
                isLoading: false,
                hasLoaded: true,
                error: null,
            });
        } catch (error) {
            const message =
                error instanceof Error ? error.message : String(error);
            set({ isLoading: false, error: message });
        }
    },
    saveScenarioItem: async (itemId, nameOverride) => {
        if (!get().hasLoaded) {
            const templates = await readCatalog();
            set({ templates: sortTemplates(templates), hasLoaded: true });
        }

        const scenarioState = useScenarioStore.getState();
        const template = createTemplateFromScenarioItem(scenarioState, itemId);
        if (!template) {
            return null;
        }
        const trimmedName = nameOverride?.trim();
        const templateToSave = trimmedName
            ? { ...template, name: trimmedName }
            : template;

        const templates = sortTemplates([...get().templates, templateToSave]);
        set({ templates, error: null });
        await writeCatalog(templates);
        return templateToSave;
    },
    addTemplates: async (incomingTemplates) => {
        if (!get().hasLoaded) {
            const templates = await readCatalog();
            set({ templates: sortTemplates(templates), hasLoaded: true });
        }

        const templatesForCatalog =
            prepareTemplatesForCatalog(incomingTemplates);
        const templates = sortTemplates([
            ...get().templates,
            ...templatesForCatalog,
        ]);
        set({ templates, error: null });
        await writeCatalog(templates);
        return templatesForCatalog.length;
    },
    updateTemplateName: async (templateId, name) => {
        if (!get().hasLoaded) {
            const templates = await readCatalog();
            set({ templates: sortTemplates(templates), hasLoaded: true });
        }

        const trimmedName = name.trim();
        if (!trimmedName) {
            return;
        }

        const updatedAt = new Date().toISOString();
        const templates = sortTemplates(
            get().templates.map((template) =>
                template.id === templateId
                    ? { ...template, name: trimmedName, updatedAt }
                    : template
            )
        );
        set({ templates, error: null });
        await writeCatalog(templates);
    },
    deleteTemplate: async (templateId) => {
        if (!get().hasLoaded) {
            const templates = await readCatalog();
            set({ templates: sortTemplates(templates), hasLoaded: true });
        }

        const templates = get().templates.filter(
            (template) => template.id !== templateId
        );
        set({ templates, error: null });
        await writeCatalog(templates);
    },
}));

export { createAssetLibraryFile, parseAssetTemplates };
