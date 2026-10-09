#!/usr/bin/env python3
"""Run Parts 2-4 against cell_counts.db and write results to outputs/."""
from pathlib import Path

import pandas as pd
import plotly.express as px

import analysis as an

OUT = Path(__file__).resolve().parent / "outputs"
OUT.mkdir(exist_ok=True)


def boxplot(freq):
    fig = px.box(
        freq, x="population", y="percentage", color="response",
        category_orders={"response": ["yes", "no"]},
        color_discrete_map={"yes": "#1f77b4", "no": "#ff7f0e"},
        points="outliers",
        labels={"percentage": "Relative frequency (%)", "population": "Cell population",
                "response": "Responder"},
        title="Melanoma, miraclib, PBMC: responders vs non-responders",
    )
    return fig


def main():
    conn = an.connect()

    summ = an.summary_table(conn)
    summ.to_csv(OUT / "part2_summary_table.csv", index=False)
    print(f"Part 2: {len(summ)} rows -> outputs/part2_summary_table.csv")

    freq = an.responder_frequencies(conn)
    res = pd.concat([an.stats_by_sample(freq), an.stats_by_subject(freq),
                     an.stats_baseline(freq)])
    res.to_csv(OUT / "part3_stats.csv", index=False)
    fig = boxplot(freq)
    fig.write_html(OUT / "part3_boxplot.html")
    try:
        fig.write_image(OUT / "part3_boxplot.png", width=1000, height=600)
    except Exception as e:  # kaleido/Chrome may be unavailable; HTML is enough
        print("PNG export skipped (kaleido not installed); HTML boxplot written.")
    print("\nPart 3 (all samples, BH-adjusted):")
    print(res[res.analysis == "all samples"][
        ["population", "median_responder", "median_nonresponder", "p_value",
         "q_value_BH", "significant_q<0.05"]].to_string(index=False))
    print("\nSensitivity analyses significant populations:")
    for a, g in res.groupby("analysis"):
        print(f"  {a}: {g.loc[g['significant_q<0.05'], 'population'].tolist()}")

    base = an.baseline_samples(conn)
    base.to_csv(OUT / "part4_baseline_samples.csv", index=False)
    print(f"\nPart 4: {len(base)} baseline melanoma/miraclib PBMC samples")
    for name, df in an.baseline_breakdowns(conn).items():
        df.to_csv(OUT / f"part4_{name}.csv", index=False)
        print(df.to_string(index=False), "\n")
    avg = an.avg_b_cell_male_responders(conn)
    (OUT / "part4_avg_b_cell.txt").write_text(f"{avg:.2f}\n")
    print(f"Avg B cells, melanoma males, responders, t=0: {avg:.2f}")


if __name__ == "__main__":
    main()
