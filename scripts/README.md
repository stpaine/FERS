# Scripts

Run these helpers from the repository root.

| Script | Purpose |
| --- | --- |
| `analyze_clang_tidy_fixes.py` | Read `fixes.yaml` and rank diagnostics by frequency; requires PyYAML. |
| `migrate_fers_xml.py` | Convert a legacy scenario to the current schema, prompting when references need disambiguation. |
| `run_tests_gen_coverage.sh` | Configure/build the coverage preset, run CTest, and generate HTML under `build/coverage/coverage_report/`. Requires `lcov` and `genhtml` for reporting. |
| `summarize_build_log.py` | Summarize compiler diagnostics and build timings, with optional Markdown or JSON output. |
| `verify_versions.py` | Check version consistency across the monorepo manifests and CMake inputs. |

```bash
python3 scripts/verify_versions.py
python3 scripts/migrate_fers_xml.py old.fersxml new.fersxml
bash scripts/run_tests_gen_coverage.sh
```

The coverage helper accepts `--tests-only` to skip report generation and `--show-passing` to display successful CTest cases.
