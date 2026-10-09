"""Analysis functions (Parts 2-4). All queries run against cell_counts.db."""
import sqlite3
from pathlib import Path

import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

DB = Path(__file__).resolve().parent / "cell_counts.db"


def connect(db=DB):
    return sqlite3.connect(db)


# ---------------- Part 2 ----------------
SUMMARY_SQL = """
SELECT cc.sample_id AS sample,
       t.total_count AS total_count,
       cc.population AS population,
       cc.count AS count,
       100.0 * cc.count / t.total_count AS percentage
FROM cell_counts cc
JOIN (SELECT sample_id, SUM(count) AS total_count
      FROM cell_counts GROUP BY sample_id) t USING (sample_id)
ORDER BY cc.sample_id, cc.population
"""


def summary_table(conn) -> pd.DataFrame:
    return pd.read_sql_query(SUMMARY_SQL, conn)


# ---------------- Part 3 ----------------
def responder_frequencies(conn) -> pd.DataFrame:
    """Relative frequencies for melanoma + miraclib + PBMC samples, with response."""
    q = """
    WITH s AS (%s)
    SELECT s.sample, sub.subject_id AS subject, sub.response,
           s.population, s.percentage, sa.time_from_treatment_start AS time
    FROM s
    JOIN samples sa  ON sa.sample_id = s.sample
    JOIN subjects sub ON sub.subject_id = sa.subject_id
    WHERE sub.condition = 'melanoma' AND sub.treatment = 'miraclib'
      AND sa.sample_type = 'PBMC' AND sub.response IN ('yes','no')
    """ % SUMMARY_SQL.replace("ORDER BY cc.sample_id, cc.population", "")
    return pd.read_sql_query(q, conn)


def _test(groups, label):
    rows = []
    for pop, g in groups:
        r = g.loc[g.response == "yes", "percentage"]
        n = g.loc[g.response == "no", "percentage"]
        u, p = stats.mannwhitneyu(r, n, alternative="two-sided")
        rows.append(dict(
            population=pop, n_responder=len(r), n_nonresponder=len(n),
            median_responder=r.median(), median_nonresponder=n.median(),
            mean_responder=r.mean(), mean_nonresponder=n.mean(),
            # rank-biserial effect size: >0 means higher in responders
            rank_biserial=2 * u / (len(r) * len(n)) - 1,
            p_value=p,
        ))
    out = pd.DataFrame(rows)
    out["q_value_BH"] = multipletests(out["p_value"], method="fdr_bh")[1]
    out["significant_q<0.05"] = out["q_value_BH"] < 0.05
    out.insert(0, "analysis", label)
    return out


def stats_by_sample(freq: pd.DataFrame) -> pd.DataFrame:
    """Mann-Whitney U per population, every sample as an observation (as specified)."""
    return _test(freq.groupby("population"), "all samples")


def stats_by_subject(freq: pd.DataFrame) -> pd.DataFrame:
    """Sensitivity analysis: average each subject's samples first, so that
    repeated measures from one patient are not treated as independent."""
    m = freq.groupby(["subject", "response", "population"], as_index=False)["percentage"].mean()
    return _test(m.groupby("population"), "subject means")


def stats_baseline(freq: pd.DataFrame) -> pd.DataFrame:
    """Sensitivity analysis: baseline (time 0) samples only -- one per subject,
    and the time point most relevant for *predicting* response."""
    return _test(freq[freq.time == 0].groupby("population"), "baseline only (t=0)")


# ---------------- Part 4 ----------------
BASELINE_SQL = """
SELECT sa.sample_id, sa.sample_type, sa.time_from_treatment_start,
       sub.subject_id, sub.project_id AS project, sub.condition, sub.treatment,
       sub.response, sub.sex, sub.age
FROM samples sa JOIN subjects sub USING (subject_id)
WHERE sub.condition = 'melanoma' AND sub.treatment = 'miraclib'
  AND sa.sample_type = 'PBMC' AND sa.time_from_treatment_start = 0
ORDER BY sa.sample_id
"""


def baseline_samples(conn) -> pd.DataFrame:
    return pd.read_sql_query(BASELINE_SQL, conn)


def baseline_breakdowns(conn) -> dict:
    b = baseline_samples(conn)
    subj = b.drop_duplicates("subject_id")
    return {
        "samples_per_project": b.groupby("project").size().rename("n_samples").reset_index(),
        "subjects_by_response": subj.groupby("response").size().rename("n_subjects").reset_index(),
        "subjects_by_sex": subj.groupby("sex").size().rename("n_subjects").reset_index(),
    }


AVG_B_SQL = """
SELECT ROUND(AVG(cc.count), 2)
FROM cell_counts cc
JOIN samples sa  ON sa.sample_id = cc.sample_id
JOIN subjects sub ON sub.subject_id = sa.subject_id
WHERE cc.population = 'b_cell' AND sub.condition = 'melanoma'
  AND sub.sex = 'M' AND sub.response = 'yes'
  AND sa.time_from_treatment_start = 0
"""


def avg_b_cell_male_responders(conn) -> float:
    return conn.execute(AVG_B_SQL).fetchone()[0]
