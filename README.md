# Teiko technical: miraclib immune-cell analysis

## Run (GitHub Codespaces)
```
make setup       # installs uv if needed, then `uv sync` (pyproject.toml)
make pipeline    # load_data.py -> cell_counts.db, then run_analysis.py -> outputs/
make dashboard   # Streamlit dashboard on port 8501
```
Dashboard link: https://upgraded-pancake-9w46jj46w7x2j5j-8501.app.github.dev/

## Layout
| File | Purpose |
|---|---|
| `load_data.py` | SQLAlchemy models: creates `cell_counts.db` and loads `cell-count.csv` |
| `analysis.py` | SQL/pandas/stats functions for Parts 2-4 |
| `run_analysis.py` | Runs the analysis, writes CSVs, boxplot and stats to `outputs/` |
| `dashboard.py` | Streamlit dashboard reading from the database |

## Database schema
- `projects(project_id)`
- `subjects(subject_id, project_id, condition, age, sex, treatment, response)`: patient-level attributes (verified constant per subject)
- `samples(sample_id, subject_id, sample_type, time_from_treatment_start)`
- `cell_counts(sample_id, population, count)`: long format, one row per sample x population

Rationale: subject attributes are stored once instead of repeated for each of a subject's
3 samples; counts in long format make the Part 2 percentages a simple group-by and let new
cell populations (or thousands of samples/projects) be added without schema changes. Indexes
cover the condition/treatment/response/sex and sample type/time filters used in Parts 3-4.

## Part 3: responders vs non-responders

Melanoma, miraclib, PBMC samples: 1,968 samples (993 responders, 975 non-responders).
Test: two-sided Mann-Whitney U per population; Benjamini-Hochberg correction across the 5 populations.

| Population | Mean % responders | Mean % non-responders | p-value | q-value (BH) |
|---|---|---|---|---|
| b_cell | 9.80 | 9.99 | 0.056 | 0.139 |
| cd8_t_cell | 24.88 | 24.94 | 0.639 | 0.639 |
| cd4_t_cell | 30.54 | 29.90 | 0.013 | 0.067 |
| nk_cell | 14.84 | 15.07 | 0.121 | 0.202 |
| monocyte | 19.94 | 20.08 | 0.163 | 0.204 |

**Conclusion:** Only `cd4_t_cell` is nominally significant (p = 0.013), with a slightly higher
share in responders. It does not survive correction for testing 5 populations (q = 0.067), and the
effect is small (about 0.6 percentage points), so I treat it as a hypothesis worth following up,
not a validated predictor. No other population differs. Sensitivity analyses (subject averages,
baseline-only) are in `outputs/part3_stats.csv`.

## Part 4 results
- 656 baseline melanoma/miraclib PBMC samples; prj1: 384, prj3: 272
- Subjects: 331 responders / 325 non-responders; 344 male / 312 female
- Average B cells, melanoma males, responders, time 0 (all sample and treatment types): **10206.15**
