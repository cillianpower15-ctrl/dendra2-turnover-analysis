# Dendra2 turnover analysis

Statistical analysis and publication figures for Dendra2 optical pulse labelling
experiments in organotypic brain slice cultures.

You give it one spreadsheet of per-slice, per-day fluorescence measurements. It
returns a time-stamped results folder containing a plain-English HTML report,
publication figures and the underlying statistical tables.

This is the second half of a two-part workflow. Image registration, z-plane
matching and extraction of the fluorescence measurements are handled by
[WittonLab/dendra2_Analysis](https://github.com/WittonLab/dendra2_Analysis)
(MATLAB), which writes the `signalData_Days.xlsx` values that go into the
workbook used here.

## What it does

- Validates the input workbook and reports the quality-control state of every
  slice before anything is fitted.
- Fits the prespecified primary model: a constrained nonlinear mixed-effects
  model of photoconverted (red) signal retention, with a slice-level random
  effect on the rate constant.
- Runs the prespecified secondary linear mixed models on the red and green
  dF/F0 trajectories.
- Produces sensitivity analyses, including leave-one-slice-out.
- Writes curated figures, tables and a report, plus a full provenance record
  (input checksum, package versions, session information).

The statistical plan is fixed in advance and is in
[`docs/Statistical_analysis_plan.md`](docs/Statistical_analysis_plan.md). Any
departure from it is a post hoc analysis and should be reported as one.

## Requirements

- R 4.1 or later. Download from <https://cran.r-project.org>.
- R packages: `readxl`, `dplyr`, `tidyr`, `ggplot2`, `lme4`, `lmerTest`,
  `purrr`, `patchwork`, `emmeans`, `nlme`. The launcher installs any that are
  missing on the first run, so be online the first time.

No other software is needed.

## Running it

1. Open `Dendra2_Input_Workbook.xlsx`. Enter one row per slice per imaging day
   on the **Master Data** sheet, check the **Analysis Settings** sheet, then
   **save and close the workbook**. Leaving it open in Excel is the most common
   cause of a failed run.
2. Start the analysis:
   - **Windows:** double-click `RUN_DENDRA2_ANALYSIS.bat` and select the
     workbook when the dialog appears.
   - **macOS or Linux:** `./run_analysis.sh Dendra2_Input_Workbook.xlsx`
   - **Any platform, from R:**
     `Rscript --vanilla scripts/run_analysis.R Dendra2_Input_Workbook.xlsx`
3. When it prints `ANALYSIS COMPLETE`, open `Dendra2_Analysis_Report.html` in
   the new results folder.

Each run writes its own time-stamped folder, so no previous run is ever
overwritten.

## Input workbook

The workbook that ships here is a blank template. The two italic rows beginning
`e.g.` are examples and are ignored by the analysis; overwrite or delete them.

**Master Data** needs seven columns as a minimum:

| Column | Meaning |
| --- | --- |
| `Group` | Experimental condition. Spelled identically on every row of that condition. |
| `Slice ID` | Repeated-measures identifier. The same slice keeps its ID across days. |
| `RecordingDay` | Imaging day. Use `1` for pre-conversion and `1.1` for the post-conversion baseline. |
| `GFP_RawSignal` | Mean raw green signal from the common registered volume. |
| `GFP_dFF0` | Green change from the structural baseline. |
| `RFP_RawSignal` | Mean raw red signal from the common registered volume. |
| `RFP_dFF0` | Red change from the structural baseline. |

Optional but strongly recommended: `PreparationID` (the animal or culture
preparation a slice came from), `Include`, `ExclusionReason`,
`RetainedZPlanes`, `RegistrationQC` and `SignalBackgroundRatio`.

Without `PreparationID` the analysis treats every slice as independent and says
so in the report. Fill it in before publishing.

`Include` is whole-slice only: a slice is either in or out, and every excluded
slice needs a reason. Day-specific exclusion within a slice is not permitted.

## Analysis Settings

Everything that changes the analysis lives on the **Analysis Settings** sheet,
not in the code. Set these before you look at any group results and do not
change them afterwards.

| Setting | Default | Notes |
| --- | --- | --- |
| Reference group | `Vehicle` | All contrasts are expressed relative to this group. |
| GFP baseline day | `1` | Defined as zero and excluded from the green model. |
| RFP baseline day | `1.1` | Structural zero; nonlinear time starts here. |
| Minimum acceptable decay R2 | `0.8` | Supporting per-slice fits below this are flagged. |
| Minimum GFP photoconversion loss fraction | `0.2` | Conversion QC. |
| Minimum RFP photoconversion fold increase | `2` | Conversion QC. |
| Minimum signal/background ratio | `2` | Measurements below this are flagged. |
| Primary analysis | Constrained nonlinear RFP rate model | Frozen endpoint. |
| Random slope strategy | Use if estimable | Intercept fallback only if singular. |
| Exclude backgrounds | *(blank)* | Optional. Drops any group whose name contains this text, comma separated. |

**`Exclude backgrounds` deserves care.** It matches on substrings and is case
insensitive, so a value of `hSyn` would drop every group whose name contains
`hSyn`, including names that use it as a promoter rather than as a condition.
Leave it blank unless you have a prespecified reason to drop a group, and check
the warning the run prints if you use it.

## Output

```
Dendra2_Results_<timestamp>/
  Dendra2_Analysis_Report.html     start here
  ANALYSIS_SUMMARY.txt             headline tests in plain text
  figures/                         curated figures, PDF and PNG, with legends
  tables/                          every statistical table as CSV
  run details/                     input snapshot, checksums, package versions, logs
```

## Reproducibility

Every run copies the input workbook, records its MD5 checksum, writes the R
version and the version of every package used, and keeps the console logs. The
results folder is therefore self-describing: it holds the inputs, the software
state and the outputs together.

## Citation

If you use this code, please cite this repository and the companion image
analysis repository.

> [[TO ADD: citation for this repository, and the paper DOI once published]]

## Licence

BSD 3-Clause. See [LICENSE](LICENSE).
