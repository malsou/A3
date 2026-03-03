# Issue Breakdown: A3 screening code

## Issue 1 — Morris elementary effects are computed on a different design than was evaluated

**Severity:** High (invalidates Morris sensitivity interpretation)

**Where:** `a3_run_screening.m`

**Problem**
`Xn_all` is generated as a Morris one-factor-at-a-time trajectory set, then denormalized to `Xraw`. After that, `enforce_opr_bounds` may randomly resample `LPR/HPR` per row. This mutates evaluated design points away from the original Morris trajectory. However, `morris_analyze` is still called with the original `Xn_all`, not with the actual evaluated points.

**Why this is wrong**
Morris elementary effects require known single-factor finite differences between consecutive trajectory points. Row-level random resampling breaks those finite-difference assumptions, so `mu*` and `sigma` are no longer valid Morris statistics.

**Evidence**
- Design is generated in normalized space and denormalized.
- OPR enforcement mutates rows.
- Analysis still uses original `Xn_all`.

**Acceptance criteria for fix**
- Either: enforce feasibility while preserving valid Morris transitions,
- Or: re-generate the analysis on the post-constraint normalized points with guaranteed one-factor transitions,
- Or: switch to a different global screening method robust to constrained resampling.

---

## Issue 2 — Mixed `VariableNames` types may break dataset table creation

**Severity:** Medium (runtime error risk)

**Where:** `a3_run_screening.m`

**Problem**
The dataset table call builds `VariableNames` using mixed types in one concatenation: `cellstr(varNames)` (cell array) and string scalars. MATLAB expects consistent container type for variable names.

**Why this is wrong**
Mixed cell/string concatenation can fail or behave inconsistently by MATLAB version.

**Acceptance criteria for fix**
- `VariableNames` provided as one consistent type (all cellstr or all string array).
- Script runs through dataset export without `array2table` naming errors.

---

## Issue 3 — Categorical comparison row assembly mixes string and numeric data in one matrix

**Severity:** Medium (data type corruption risk)

**Where:** `a3_run_screening.m`

**Problem**
`catRows` is grown by concatenating categorical strings (`nozzle`, `cooling`) with numeric values in a single horizontal vector. In MATLAB, this tends to coerce to a common type (often string), converting numeric values to text.

**Why this is wrong**
The resulting table can lose numeric column types, breaking numeric post-processing and statistics.

**Acceptance criteria for fix**
- Build `Tcat` with typed columns (e.g., table columns assigned separately or via struct arrays).
- Numeric outputs remain numeric in `categorical_compare.csv`.

---

## Issue 4 — OPR enforcement may silently leave infeasible rows

**Severity:** Low/Medium (silent quality issue)

**Where:** `enforce_opr_bounds.m`

**Problem**
After up to 200 retries, the function exits the loop without warning/error even if a row still violates OPR bounds.

**Why this is wrong**
Silent failure can contaminate analysis while appearing successful.

**Acceptance criteria for fix**
- Emit warning/error when retry limit is reached and row remains infeasible,
- And/or guarantee feasibility by deterministic projection/repair.

