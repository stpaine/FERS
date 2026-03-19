// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

/**
 * A utility function for exhaustive checks in switch statements.
 *
 * If a switch statement should handle all possible cases of a union type,
 * you can use this function in the `default` case. TypeScript will ensure
 * that the type of the variable passed to this function is `never`, meaning
 * all other cases have been exhausted. If you add a new type to the union
 * and forget to handle it in the switch, the compiler will raise an error
 * at the call site of this function.
 *
 * @param x - The variable that should be of type `never`.
 * @returns This function never returns; it throws an error.
 */
export function assertNever(x: never): never {
    throw new Error(`Unexpected object: ${JSON.stringify(x)}`);
}

/**
 * Creates a new object by omitting specified keys from an existing object.
 * This is a type-safe alternative to using destructuring with rest syntax
 * when trying to avoid `no-unused-vars` linting errors.
 *
 * @param obj The source object.
 * @param keys The keys to omit from the new object.
 * @returns A new object without the specified keys.
 */
export function omit<T extends Record<string, unknown>, K extends keyof T>(
    obj: T,
    ...keys: K[]
): Omit<T, K> {
    const ret = { ...obj };
    for (const key of keys) {
        delete ret[key];
    }
    return ret;
}
