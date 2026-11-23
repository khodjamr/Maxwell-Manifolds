# Maxwell tuning manifolds reproducibility package

Code, data, and verification materials accompanying **Maxwell Tuning
Manifolds: Quadratic Structure and Critical Transitions**, by Mohamed R.
Khodja.

The package contains the data for all ten manuscript figures, the native
Mathematica 9 numerical package and launchers, publication figure PDFs, exact
rational certificate inputs, and independent verification programs. It also
retains the auxiliary tables that support statements in the manuscript.

## Choose the task you want to reproduce

| Task | Entry point | Requirements | Output |
|---|---|---|---|
| Replot the supplied data in the manuscript layout | `Maxwell_Manifold_clean_figures.nb` | Mathematica 9 | Ten displayed plots and PDFs in `clean_figures/` |
| Recompute the numerical Maxwell figure data | `RunAllFigures.nb` | Mathematica 9; the adjacent `.m` package | Eleven primary CSVs in `data/`, ten original-style PDFs in `figures/` |
| Check download integrity, data summaries and all exact certificates | `python analysis/verify_release.py` | Python 3.8 or later, standard library only | Read-only PASS/FAIL report |
| Audit the supplied numerical values independently | `python analysis/audit_numerics.py` | Optional scientific packages below | `verification_results/numerical_consistency.json` |
| Run the independent numerical reference implementation | `python analysis/reproduce.py` | Optional scientific packages below | A separate `independent_run/` tree |

The Mathematica numerical package is **3.2.5**. The CSV-only plotting source is
**1.0.11**. These are different programs with separate version numbers.

## Replot without recomputing Maxwell fields

1. Extract the entire folder. Keep the notebook beside the `data` directory.
2. Open `Maxwell_Manifold_clean_figures.nb` in Mathematica 9.
3. Start a fresh local kernel and choose **Evaluation > Evaluate Notebook**.
   Allow the notebook's initialization cell to evaluate if prompted.
4. The notebook reads eleven CSV tables, displays all ten figures and exports
   the PDFs to `clean_figures/`. Each panel's legend-width control remains
   beside its figure command.

This notebook does not load the numerical `.m` package, solve a Maxwell
problem or verify an interval certificate. Its `mcfSelfTest[]` checks the CSV
inputs and plotting requirements only. The publication PDFs are also
provided in [clean_figures/](clean_figures/).

## Recompute the Maxwell figure data

1. Keep `RunAllFigures.nb` and `MaxwellManifoldFeasibility.m` in the same folder.
2. Start a fresh Mathematica 9 local kernel and open `RunAllFigures.nb`.
3. Evaluate the loading cell and `SelfTest[]`. Continue after the test returns
   `True` and prints `SelfTest PASS v3.2.5`.
4. Evaluate the figure sections in order, or evaluate the whole notebook.
   The spheroidal calculation and critical-point scan can take substantial
   time. Returned plots are displayed in the launcher.
5. To obtain the manuscript styling from the regenerated data, evaluate the
   CSV-only notebook afterwards.

The numerical package creates `figures/` and `data/` when needed. Re-running a
figure replaces its corresponding output files. Individual commands are
listed in `FIGURE_SCRIPT_MAP.md`; for example, `RunFigure7[]` and
`RunWeakLossStudy[]` generate Figures 7 and 10 separately after a successful
load and self-test. Neither program depends on cached notebook output.

The primary figure data were generated with native Mathematica 9.0 on
Windows. A native diagnostic records a passed self-test for numerical package
3.2.5. The independent verification did not include a fresh native Mathematica
run. Execution provenance, verification results, and package-history details
are documented in [AUDIT_SUMMARY.md](AUDIT_SUMMARY.md) and
[audit/software_environment.json](audit/software_environment.json).

## Exact verification

From the extracted folder run:

```text
python analysis/verify_release.py
```

This read-only command verifies the manifest, primary and auxiliary CSV
counts, reported fits, candidate counts, fixed-power resource range, two
critical-point certificates, quantitative dipole bounds, and loss-gradient
signs. It also encloses the two magnetic composition-shift derivatives used
to normalize Figure 10(b). No NumPy, SciPy, mpmath, SymPy or Mathematica is
needed for this command. Do not run proof verification with Python's `-O` or
`-OO` flags, which disable assertions; the interval engine rejects that mode.

The individual exact replayers are:

```text
python -S analysis/verify_certificates.py
python -S analysis/small_size_certificate.py
python -S analysis/verify_loss_direction.py
```

The `-S` flag demonstrates that no installed scientific packages are needed.
Unlike the combined read-only checker, these individual commands rewrite
their corresponding readable result files in `audit/`.

### What the certificate file contains

`audit/critical_point_certificates.json` is **one input file containing two
certificate records**, one saddle and one definite example. The exact rational
centres, box radii and preconditioners are essential proof inputs.
Verification recomputes the interval bounds instead of
trusting stored Boolean flags or decimal summaries.

`analysis/rational_interval.py` supplies shared fraction-based interval and
Taylor-jet arithmetic. Its 360-bit dyadic rounding and 91-term elementary
function bounds match Appendix E. Floating-point candidate selection is
separate, in `analysis/rational_certificate.py`; generating candidates is not
needed to replay the supplied proofs. That optional generator writes to
`verification_results/generated_certificates/` and does not replace the
archived certificate inputs.

The critical-point checks certify local existence and uniqueness inside the
two supplied boxes. They do not certify every candidate in Figure 8. The
small-size checks verify sufficient bounds on the stated domains, not a sharp
critical-size threshold. The loss check certifies the nonzero full-material
gradient at the two lossless points; finite-loss coordinates remain numerical
continuations.

Decimal summaries are for reading. Rational endpoints drive interval
inclusion and contraction. Rounded Table II coordinates are not the exact
centres of the proof boxes.

## Independent numerical and symbolic checks

The optional reference implementation uses Python scientific packages. Install
the audit versions, if needed, with:

```text
python -m pip install -r requirements-independent.txt
```

Use Python 3.11 or later for these pinned scientific packages. The lighter
standard-library proof replay above needs Python 3.8 or later.

To check the supplied CSV values without replacing them:

```text
python analysis/audit_numerics.py
```

To run exact and symbolic checks in an isolated reference tree:

```text
python analysis/reproduce.py --checks-only
```

To regenerate the full independent study, its auxiliary tables and its own
figures:

```text
python analysis/reproduce.py
```

An alternative output directory can be selected with `--output PATH`. The
driver copies the code and required inputs there before execution and checks
that canonical Mathematica data and publication PDFs remain unchanged.

The independent implementation uses some different grids and CSV column
names. Its output should not replace `data/` as input to the publication
plotting notebook. Agreement of physical quantities, convergence and
classification is the relevant comparison; byte-identical plots are not
expected. The independent full workflow is an additional validation route,
not a prerequisite for plotting or replaying the exact certificates.

## Data and evidence map

| Location | Role |
|---|---|
| `data/fig1_...csv` through `data/fig10_...csv` | Eleven primary plotting tables; Figure 8 has two tables |
| `data/fig4_extended_asymptotics.csv` | Five additional asymptotic cases discussed in Section V |
| `data/fig8_candidate_counts.csv` | Counts, including zero counts, on the 21-size non-exhaustive scan |
| `data/higher_multipole_persistence.csv` | Eight quadrupole perturbation tests discussed in Section VI |
| `data/high_precision_critical_points.json` | High-precision numerical working points for the independent calculations |
| `audit/critical_point_certificates.json` | Unchanged exact proof inputs for both critical points |
| `audit/certificate_replay.json` | Replayed critical-point bounds |
| `audit/small_size_certificate.json` | Replayed sufficient dipole bounds, including a failed larger box |
| `audit/loss_direction_certificate.json` | Loss-gradient and magnetic control-rate enclosures |
| `audit/numerical_consistency.json` | Checks of the supplied native figure data |
| `audit/independent_maxwell_audit.json` | Independent boundary-matching and volume-work checks |
| `audit/spheroid_check.json` | Independent double-precision spheroidal and quadrature comparisons |
| `audit/weak_loss.json` | Independently recomputed weak-loss continuations, derivatives and Hessians |
| `audit/modified_bessel_audit.json` | Negative-internal-parameter determinant checks in Appendix A |
| `audit/asymptotic_extended_cases.json` | Additional asymptotic fitted slopes |
| `audit/symbolic_series_check.json` | Exact dipole series and power-coefficient checks |
| `audit/geometry_structure_audit.json` | Symbolic operator/lift consistency checks; the analytic proof is in the paper |
| `audit/software_environment.json` | Execution provenance and dependency versions |
| `clean_figures/` | Ten publication PDFs, with drawing streams matching the manuscript |

Empty entries in the multi-panel Figure 6 CSV mark quantities that do not
belong to that dataset; they are not missing measurements. Row samples and
branch identifiers separate curves. The stored residuals were evaluated
before export at higher precision than the 17-digit CSV coordinates. Reusing
the rounded coordinates can therefore give larger residuals than the stored
ones; this does not replace or invalidate the exact certificate checks.

The Figure 4 fits use the first seven size samples. Figure 7 fits use the
first six control increments for each example. The additional asymptotic
fits use the first six samples of each case. Those fitting windows matter
when reproducing the numerical exponents quoted in the paper.

The extremely small held-out power discrepancy in the native spheroidal
data is an empirical algebraic consistency check. A separate double-precision
solver reaches a different rounding floor. Neither residual is a continuum
error certificate.

## Integrity

`MANIFEST_SHA256.csv` records every distributed file except itself. Run the
read-only release checker before regenerating outputs to test the downloaded
snapshot. New PDF metadata and floating-point output can legitimately change
checksums after regeneration; checksums establish snapshot integrity, not
bitwise reproducibility across different software installations.

## Citation

Please cite the accompanying manuscript and the specific release of this
reproducibility package used in your work. Citation metadata are provided in
[CITATION.cff](CITATION.cff). The release identifier and checksum manifest
identify the computational materials used for reproduction.
