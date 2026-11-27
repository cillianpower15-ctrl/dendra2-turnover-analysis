# Statistical analysis plan

Pre-specified before the analysis was run. Fixed thereafter; any departure from
it is a post hoc analysis and should be reported as one.

## Experimental unit and hierarchy

`Slice ID` is the repeated-measures unit. `PreparationID` records the animal or
preparation from which a slice arose. Slice and preparation counts will be
reported separately. When complete and crossed with group, PreparationID is a
fixed blocking factor. If preparation is missing or confounded with group, the
software will flag that the group effect cannot be cleanly separated from
preparation.

## Inclusion and exclusion

Inclusion is assigned to the whole slice before analysis. A reason is mandatory
for every excluded slice. The package does not permit day-specific inclusion
within a slice. Excluded rows are retained in `00_excluded_records.csv`.

## Primary endpoint and test

The primary endpoint is the post-photoconversion RFP turnover-rate constant
`K`. Fraction retained is fitted as `(1-P)exp(-Kt)+P`, with `K=exp(logK)` and
`P=plogis(logitP)`. Slice receives a random effect on `logK`.

The structural post-conversion baseline defines fraction retained as one at
`t=0`, but it is excluded from the likelihood. Model M1 has a shared rate; M2
adds group effects on `logK`. The prespecified primary test is the ML
likelihood-ratio comparison M1 versus M2 at two-sided alpha 0.05. M3 adds group
effects on the plateau and is secondary. Pairwise rate ratios are Holm
adjusted.

## Secondary analyses

RFP and GFP dF/F0 trajectories are assessed with linear mixed models containing
Day, Group and Day-by-Group. A slice-level random slope for Day is
prespecified and retained whenever convergent and non-singular; a random
intercept is used only as an estimability fallback. Estimated trend and level
contrasts use Tukey adjustment.

Supporting slice-level summaries include constrained plateau half-life,
no-plateau and log-linear alternatives, and trapezoidal AUC from each
channel-specific structural baseline. Parametric and rank tests, effect sizes,
model sensitivity and leave-one-slice-out results are reported.

## QC and reporting

Photoconversion, retained Z planes, registration, signal/background, missing
timepoints and exclusions are reported without post hoc deletion. A non-PASS
readiness item must be resolved or justified. Estimates, 95% confidence
intervals, exact adjusted p-values, slice/preparation counts and all exclusions
will be reported. Analysis software versions, input/output checksums and full
session information are archived.
