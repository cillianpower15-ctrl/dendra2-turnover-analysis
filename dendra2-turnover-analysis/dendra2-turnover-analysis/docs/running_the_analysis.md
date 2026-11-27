# Running the analysis

How to run the pipeline, what each script does, what comes out, and what to do
when something fails.

## Before the first run

- Install R from <https://cran.r-project.org>. The default options are fine.
- Be connected to the internet for the first run. Any missing R packages
  install themselves. After that it works offline.

## Each run

1. Open `Dendra2_Input_Workbook.xlsx`, enter or paste your measurements on the
   **Master Data** sheet, confirm the **Analysis Settings** sheet, then **save
   and close the workbook**. Leaving it open in Excel is the most common reason
   a run fails.
2. Start the analysis:
   - Windows: double-click `RUN_DENDRA2_ANALYSIS.bat` and pick the workbook.
   - macOS or Linux: `./run_analysis.sh Dendra2_Input_Workbook.xlsx`
   - Any platform: `Rscript --vanilla scripts/run_analysis.R <workbook>`
3. Leave the console open. It prints each of the eight steps and takes a few
   minutes. When it prints `ANALYSIS COMPLETE`, open
   `Dendra2_Analysis_Report.html` in the new results folder.

Each run writes its own time-stamped folder, so nothing is ever overwritten.

## Entering data

One row per slice per imaging day. `Group` names must be spelled identically on
every row of the same condition; two spellings become two groups.

Use `RecordingDay` `1` for the pre-conversion image and `1.1` for the
post-conversion image acquired immediately afterwards. Later days are the
elapsed days of the chase.

Adding groups needs no code change. The pipeline takes the group list from the
sheet. If the names follow a `Background <something> Construct` pattern, the
by-background figures and within-background statistics appear automatically.

## Choosing which groups are analysed

Two independent mechanisms exclude data, and they do different jobs.

- **`Include` column (per slice).** Set to `No` to drop one slice, and give a
  reason in `ExclusionReason`. A reason is mandatory. Inclusion is whole-slice:
  you cannot drop a single day from a slice.
- **`Exclude backgrounds` setting (per group).** Drops any group whose name
  *contains* the given text, comma separated, case insensitive. Leave it blank
  to analyse every group.

`Exclude backgrounds` matches substrings, so it will also drop groups whose
names merely contain that text for another reason, such as a promoter name.
Leave it blank unless you have a prespecified reason, and read the warning the
run prints when it fires.

## Inside `scripts`

| File | Role |
| --- | --- |
| `run_analysis.R` | Orchestrator. Runs the numbered scripts in order. |
| `00_setup.R` | Shared configuration, data loading and validation, figure theme. Reads settings from the workbook. |
| `01_statistics.R` | Quality control, linear mixed models, per-slice metrics. |
| `02_figures.R` | Decay figures. |
| `03_sensitivity_checks.R` | Robustness checks, including leave-one-slice-out. |
| `04_nonlinear_mixed_model.R` | The primary analysis: the turnover-rate model. |
| `06_publication_figure.R` | Composite main figure. |
| `07_factorial_figure.R` | By-background figure with per-day p-values. |
| `08_curate_figures.R` | Trims and renames figures down to the publication set. |
| `05_build_report.R` | Writes the HTML report last, once every figure exists. |

Steps run in the order `01, 02, 03, 04, 06, 07, 08, 05`, so the report is
written last and can reference every figure that exists on disk.

## What comes out

| Output | What it is |
| --- | --- |
| `Dendra2_Analysis_Report.html` | Read this first. Everything on one page. |
| `ANALYSIS_SUMMARY.txt` | The headline numbers in plain text. |
| `figures/Figure1_MAIN` | Decay by background with per-day p-values, and the rate ratios. Panels also saved separately. |
| `figures/FigureS1_per_replicate_metrics` | Half-life and AUC per slice. |
| `figures/FigureS2_GFP_synthesis` | Unconverted-channel recovery by background. |
| `figures/Figure1_MAIN_LEGEND.txt` | Full legend with the current numbers, plus a shorter plain-language version. |
| `tables/` | The exact numbers behind every figure, as CSV. |
| `run details/` | Input snapshot and checksum, package versions, methods scaffold, full statistical transcripts. |

A fresh run regenerates all of this in a new folder, so anything deleted by
accident comes back by running again.

## If something goes wrong

| Symptom | Fix |
| --- | --- |
| "Could not read the data", or a permission error | The workbook is still open in Excel. Close it and run again. |
| A message about missing R packages | Connect to the internet and run again. They install themselves. |
| "No usable rows were found in 'Master Data'" | The sheet is empty, or every row still starts with the `e.g.` example text. |
| "Reference group not in this sheet" | Set `Reference group` on Analysis Settings to your exact control name. |
| A group you expected is missing | Check `Exclude backgrounds` on Analysis Settings, and the `Include` column for those slices. |
| One condition appears as two groups | Two spellings of the same `Group` name. Make them identical. |
| "Every excluded slice needs a prespecified ExclusionReason" | Fill in `ExclusionReason` for every row where `Include` is `No`. |
| "Include must be identical across every day of a slice" | A slice is marked `Yes` on some days and `No` on others. Pick one. |

## Reporting the analysis

A methods scaffold is written to `run details/METHODS_REPORTING_TEMPLATE.txt`
on every run, pre-filled with the software versions and settings actually used.
Start from that rather than from memory, and read it alongside
[`Statistical_analysis_plan.md`](Statistical_analysis_plan.md).
