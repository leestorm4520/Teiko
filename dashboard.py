"""Streamlit dashboard for Parts 2-4. Run with `make dashboard`."""
import pandas as pd
import plotly.express as px
import streamlit as st

import analysis as an

st.set_page_config(page_title="Loblaw Bio - miraclib trial", layout="wide")

if not an.DB.exists():
    st.error("cell_counts.db not found. Run `make pipeline` first.")
    st.stop()


@st.cache_data
def load():
    conn = an.connect()
    summ = an.summary_table(conn)
    freq = an.responder_frequencies(conn)
    stats = pd.concat([an.stats_by_sample(freq), an.stats_by_subject(freq),
                       an.stats_baseline(freq)])
    base = an.baseline_samples(conn)
    breakdowns = an.baseline_breakdowns(conn)
    avg = an.avg_b_cell_male_responders(conn)
    conn.close()
    return summ, freq, stats, base, breakdowns, avg


summ, freq, stats, base, breakdowns, avg = load()

st.title("Loblaw Bio: miraclib immune-cell analysis")
t2, t3, t4 = st.tabs(["Part 2: Frequencies", "Part 3: Responders vs non-responders",
                      "Part 4: Baseline subset"])

with t2:
    st.subheader("Relative frequency of each cell population per sample")
    pops = st.multiselect("Population", sorted(summ.population.unique()),
                          default=sorted(summ.population.unique()))
    q = st.text_input("Filter by sample id (substring)")
    view = summ[summ.population.isin(pops)]
    if q:
        view = view[view["sample"].str.contains(q, case=False)]
    st.caption(f"{len(view):,} rows")
    # column_config instead of Styler: Styler caps rendering at 262,144 cells
    st.dataframe(view, width="stretch", height=500, hide_index=True,
                 column_config={"percentage": st.column_config.NumberColumn(format="%.2f")})

with t3:
    st.subheader("Melanoma, miraclib, PBMC samples")
    c1, c2 = st.columns(2)
    times = c1.multiselect("Time from treatment start (days)", [0, 7, 14], default=[0, 7, 14])
    f = freq[freq.time.isin(times)]
    fig = px.box(f, x="population", y="percentage", color="response",
                 category_orders={"response": ["yes", "no"]},
                 color_discrete_map={"yes": "#1f77b4", "no": "#ff7f0e"},
                 labels={"percentage": "Relative frequency (%)", "population": "Cell population",
                         "response": "Responder"})
    st.plotly_chart(fig, width="stretch")

    st.markdown("**Statistics** (two-sided Mann-Whitney U per population; "
                "Benjamini-Hochberg correction across the 5 populations)")
    cols = ["population", "n_responder", "n_nonresponder", "median_responder",
            "median_nonresponder", "rank_biserial", "p_value", "q_value_BH",
            "significant_q<0.05"]
    for name, g in stats.groupby("analysis", sort=False):
        st.markdown(f"*{name}*")
        st.dataframe(g[cols].style.format(
            {c: "{:.4g}" for c in ["median_responder", "median_nonresponder",
                                   "rank_biserial", "p_value", "q_value_BH"]}),
            width="stretch", hide_index=True)
    sig = stats[(stats.analysis == "all samples") & stats["significant_q<0.05"]].population.tolist()
    nominal = stats[(stats.analysis == "all samples") & (stats.p_value < 0.05)].population.tolist()
    st.info(f"Significant after multiple-testing correction (q < 0.05): "
            f"{', '.join(sig) if sig else 'none'}. "
            f"Nominal p < 0.05 before correction: {', '.join(nominal) if nominal else 'none'}.")

with t4:
    st.subheader("Baseline (time 0) melanoma PBMC samples from miraclib-treated patients")
    st.metric("Matching samples", len(base))
    c1, c2, c3 = st.columns(3)
    c1.markdown("**Samples per project**")
    c1.dataframe(breakdowns["samples_per_project"], hide_index=True)
    c2.markdown("**Subjects by response**")
    c2.dataframe(breakdowns["subjects_by_response"], hide_index=True)
    c3.markdown("**Subjects by sex**")
    c3.dataframe(breakdowns["subjects_by_sex"], hide_index=True)
    st.metric("Avg B cells: melanoma males, responders, time 0 (all sample/treatment types)",
              f"{avg:.2f}")
    with st.expander("Matching samples"):
        st.dataframe(base, width="stretch")
