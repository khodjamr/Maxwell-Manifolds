# Figure data and generator map

The eleven primary CSVs below are the supplied native Mathematica data.
The canonical publication PDFs are in `clean_figures/`. The numerical `.m`
package generates its own original-style PDFs in `figures/`; the CSV-only
notebook supplies the publication layout afterwards.

| Figure | Primary CSV in `data/` | Data rows | Native command | Optional independent source |
|---|---|---:|---|---|
| 1 | `fig1_spheroid_convergence.csv` | 11 | `RunSpheroidStudy[]` | `analysis/spheroid_check.py` |
| 2 | `fig2_static_profiles.csv` | 1266 | `RunSphericalStudy[]` | `analysis/build_study.py`, `figure2()` |
| 3 | `fig3_regular_tuning_graph.csv` | 241 | `RunSphericalStudy[]` | `analysis/build_study.py`, `figure3()` |
| 4 | `fig4_asymptotic_convergence.csv` | 15 | `RunSphericalStudy[]` | `analysis/build_study.py`, `figure4()` |
| 5 | `fig5_saddle_reconnection.csv` | 9518 | `RunSphericalStudy[]` | `analysis/build_study.py`, `figure5()` |
| 6 | `fig6_closed_loop_birth.csv` | 5519 | `RunSphericalStudy[]` | `analysis/build_study.py`, `figure6()` |
| 7 | `fig7_critical_scaling.csv` | 26 | `RunFigure7[]`, also included in `RunSphericalStudy[]` | `analysis/build_study.py`, `figure7()` |
| 8(a) | `fig8_discriminant_root_curves.csv` | 17407 | `RunCriticalMap[]` | `analysis/critical_map.py` |
| 8(b) | `fig8_critical_map.csv` | 201 | `RunCriticalMap[]` | `analysis/critical_map.py` |
| 9 | `fig9_design_interpretation.csv` | 644 | `RunDesignStudy[]` | `analysis/extended_study.py`, `figure9()` |
| 10 | `fig10_weak_loss.csv` | 12 | `RunWeakLossStudy[]` | `analysis/weak_loss.py`; independent plots/checks in `plot_weak_loss.py` |

All native commands require the adjacent package to be loaded and its
`SelfTest[]` to pass. `RunAllFigures.nb` runs these commands and displays their
plots. `Maxwell_Manifold_clean_figures.nb` reads the listed CSVs and displays
and exports all ten publication figures without repeating the calculations.

The following auxiliary tables are retained because they support manuscript
statements, even though the primary native launcher does not export them:

| Auxiliary table | Rows | Reproduction source | Manuscript role |
|---|---:|---|---|
| `fig4_extended_asymptotics.csv` | 55 | `extended_study.py`, `broadened_asymptotics()` | Additional cases and fitted slope range |
| `fig8_candidate_counts.csv` | 21 | Counts reconstructed from the primary map on its prescribed scan grid; also produced by `critical_map.py` | Counts and the empty sampled size range |
| `higher_multipole_persistence.csv` | 8 | `higher_mode_persistence.py` | Quadrupole perturbation shifts and fourth-order residual trend |

Use the isolated `analysis/reproduce.py` driver for the independent sources.
Their alternative grids and some column names differ from the canonical
Mathematica output; they are not replacement input files for the clean-figure
notebook. See `README.md` for the distinction between plotting, numerical
validation and exact interval proof.
