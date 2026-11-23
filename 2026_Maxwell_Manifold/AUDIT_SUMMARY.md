# Reproducibility audit summary

Audit date: 23 November 2026. Numerical package: 3.2.5. CSV-only plotting
source: 1.0.11. The supplied 47-page manuscript was the comparison document.

## Findings

The checked numerical results support the manuscript's reported values and
classifications. The two exact critical-point proofs, the quantitative dipole
exclusion bounds and the loss-gradient signs were replayed successfully with
standard-library exact arithmetic. The full independent reference workflow
also completed successfully in an isolated output tree.

The native Mathematica numerical package, numerical launcher, eleven primary
CSV files, three auxiliary CSV files and high-precision working-point inputs
were preserved unchanged. The two rational certificate records were preserved
unchanged in their single input file.

The publication PDFs for Figures 1 and 2(a) had legend-width differences from
the manuscript. Those two layout settings were aligned with the manuscript;
no curves, values, labels or scientific calculations changed. The ten
publication PDF drawing streams now match their manuscript figures exactly.
Saved output and local machine paths were removed from the CSV-only notebook.

## Representative reproduced quantities

| Quantity | Reproduced result |
|---|---|
| Saddle contraction bound | 3.0409905342e-14 |
| Definite contraction bound | 1.0400026932e-14 |
| Joint-domain derivative lower bound, rho <= 0.01, abs(t) <= 1/4, abs(mu) <= 100 | 0.1132962433, greater than the stated 0.113 |
| Joint-domain derivative lower bound, rho <= 0.1, abs(t) <= 1/4, abs(mu) <= 10 | 0.0092385123, greater than the stated 0.009 |
| Original bounded-contrast derivative lower bound, rho <= 0.005 | 0.1980632506, greater than the stated 0.198 |
| Figure 3 sampled scaled permittivity derivative | 0.25258496 to 0.26129649 |
| Figure 4 power-error slopes, first seven samples | 1.99728469 and 4.00020454 |
| Figure 4 gradient-error slopes, first seven samples | 2.00019854 and 4.00022697 |
| Figure 7 saddle displacement/inverse-sensitivity slopes, first six samples | 0.49995151 and -0.50001359 |
| Figure 7 definite displacement/inverse-sensitivity slopes, first six samples | 0.50000153 and -0.49999739 |
| Figure 8 candidates at rho = 3 | 21 |
| First detected size on the recorded, non-exhaustive grid | rho = 1.125 |
| Figure 9 larger-loop fixed-power resource ratio | 0.9957821013 to 1.0046762087 |
| Figure 9 larger-loop relative resource variation | -0.42178987% to +0.46762087% |
| Figure 10 zero-loss magnetic composition derivatives | 16.30703156275 and -0.02851444687408 |

All twelve finite-loss rows were independently checked using conjugate-profile
currents, boundary matching, volume work and material derivatives. The
quadrupole and additional asymptotic tables were independently regenerated
and agree with the supplied tables at their recorded precision.

## Scope of the checks

`audit/numerical_consistency.json` records which checks cover all rows and
which use samples. In particular, the high-precision root-curve, normal-form
and resource comparisons use explicit recorded subsets; all primary tables
were checked for their structure and relevant basic invariants. The exact
certificates cover only their two rational boxes. Non-exhaustive scans and
empirical spheroidal convergence are not promoted to existence or continuum
error certificates.

Native Mathematica 9 was not installed in the audit environment. The supplied
native diagnostic recorded a passed version-3.2.5 self-test; the current
audit independently checked the supplied data and ran the separate reference
workflow. No fresh native-kernel execution is claimed.

The old checksum table and version-3.0 summaries were obsolete for the current
native data. They were replaced by the current manifest and verification
records. Historical compilation logs, bibliographic lookup files, exploratory
candidate lists, the old original-style PDF/PNG set and the Figure-10-only
repair instructions were removed from this publication-oriented snapshot.

The additional asymptotic, candidate-count and quadrupole CSV files are not
obsolete: they support statements still present in the manuscript. Their
generation routes are documented in `FIGURE_SCRIPT_MAP.md`.
