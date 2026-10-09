# Teiko technical: miraclib immune-cell analysis

## Run (GitHub Codespaces)
```
make setup       # install dependencies (requirements.txt)
make pipeline    # load_data.py -> cell_counts.db, then run_analysis.py -> outputs/
make dashboard   # Streamlit dashboard on port 8501
```
Dashboard link: `http://localhost:8501` (in Codespaces, open the forwarded port 8501).

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

## Part 3 statistical approach
Population: melanoma, miraclib, PBMC samples with a known response. Two-sided
Mann-Whitney U test per population (relative frequencies are bounded and not normally
distributed), Benjamini-Hochberg FDR across the 5 populations, plus a rank-biserial effect size.
Because each subject contributes up to 3 samples, I also report two sensitivity analyses:
subject-averaged frequencies, and baseline (t=0) only.

**Result:** no population differs significantly between responders and non-responders after
multiple-testing correction in any analysis. `cd4_t_cell` has a nominal p ≈ 0.013 (q ≈ 0.067)
with a small effect (median 30.2% vs 29.7%); this should be treated as a hypothesis to
follow up, not a finding that predicts response.

## Part 4 results
- 656 baseline melanoma/miraclib PBMC samples; prj1: 384, prj3: 272
- Subjects: 331 responders / 325 non-responders; 344 male / 312 female
- Average B cells, melanoma males, responders, time 0 (all sample and treatment types): **10206.15**
