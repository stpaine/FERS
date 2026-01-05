# Security Policy

The FERS project takes security seriously. As a user, you should feel confident that the software is safe to use. This
document outlines how to report security vulnerabilities.

## Supported Versions

As this project is in active development, there are no official versioned releases yet. We provide security support for
the most recent commit on the `main` branch only.

| Branch | Supported          |
| ------ | ------------------ |
| `main` | :white_check_mark: |
| Other  | :x:                |

## Reporting a Vulnerability

> [!WARNING]
> **Do not report security vulnerabilities through public GitHub issues, pull requests, or discussions.**

If you believe you have found a security vulnerability, please report it privately through **GitHub's private
vulnerability reporting** feature. This ensures the issue can be investigated and addressed before it becomes public
knowledge.

**➡️ [Report a vulnerability here](https://github.com/davidbits/FERS/security/advisories/new)**

You can find detailed instructions on the process in GitHub's
documentation: "[Privately reporting a security vulnerability](https://docs.github.com/en/code-security/security-advisories/guidance-on-reporting-and-writing/privately-reporting-a-security-vulnerability)".

### What to Include in Your Report

To help resolve the issue quickly, please include:

- A clear description of the vulnerability and its potential impact.
- Step-by-step instructions to reproduce the issue.
- A proof-of-concept, such as a sample `.fersxml` file, that demonstrates the vulnerability.

## Scope of this Policy

Due to the nature of FERS as an offline simulation tool, our security concerns are focused on the safe handling of
user-provided input files.

### In Scope (Considered a Security Vulnerability)

- **Arbitrary Code Execution (ACE)** when parsing a malicious `.fersxml` or other input file.
- **Denial of Service (DoS)**, such as a crash or significant resource exhaustion (CPU/memory), caused by a specially
  crafted input file.
- **Path Traversal** vulnerabilities that allow reading or writing files outside of the intended scope.
- Vulnerabilities in the `fers-ui` application that could lead to local system compromise.

### Out of Scope (Considered a Bug, Not a Security Vulnerability)

- The simulator crashing on a malformed or invalid (but not maliciously crafted) `.fersxml` file. Please report these as
  regular [bug reports](https://github.com/davidbits/FERS/issues/new/choose).
- The simulator producing scientifically incorrect results. This is a correctness issue and should also be reported as a
  bug.

## The Process & Our Commitment

1. After you submit a report, it will be triaged and investigated.
2. We will keep you informed of our progress through the private advisory.
3. Upon disclosure, we will publicly credit you for your discovery unless you prefer to remain anonymous.

**Please note:** This is a solo-maintained project. We aim to acknowledge receipt of your report within **one week**,
but response times may vary. Your patience is greatly appreciated.

We do not currently offer a bug bounty program.
