import globals from 'globals';
import pluginJs from '@eslint/js';
import tseslint from 'typescript-eslint';
import pluginReact from 'eslint-plugin-react';
import pluginReactHooks from 'eslint-plugin-react-hooks';
import pluginReactRefresh from 'eslint-plugin-react-refresh';
import eslintConfigPrettier from 'eslint-config-prettier';

export default [
    // Global ignore patterns
    {
        ignores: ['dist/', 'src-tauri/', 'node_modules/', '.vite/'],
    },

    // Base configurations for all JS/TS files
    pluginJs.configs.recommended,
    ...tseslint.configs.recommended,

    // Configuration for React files
    {
        files: ['src/**/*.{js,jsx,ts,tsx}'],
        plugins: {
            react: pluginReact,
            'react-hooks': pluginReactHooks,
            'react-refresh': pluginReactRefresh,
        },
        languageOptions: {
            parserOptions: {
                ecmaFeatures: { jsx: true },
            },
            globals: {
                ...globals.browser,
            },
        },
        settings: {
            react: {
                version: 'detect', // Automatically detect the React version
            },
        },
        rules: {
            ...pluginReact.configs.recommended.rules,
            ...pluginReactHooks.configs.recommended.rules,
            'react/react-in-jsx-scope': 'off',
            'react/jsx-uses-react': 'off',
            'react/no-unknown-property': 'off',
            'react-refresh/only-export-components': [
                'warn',
                { allowConstantExport: true },
            ],
        },
    },

    // Configuration for project config files
    {
        files: ['vite.config.ts', 'eslint.config.js', 'prettier.config.js'],
        languageOptions: {
            globals: {
                ...globals.node,
            },
        },
    },

    // Prettier config must be last to override other styling rules
    eslintConfigPrettier,
];
