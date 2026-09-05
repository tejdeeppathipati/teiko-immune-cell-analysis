"""Interactive Streamlit dashboard for the immune-cell analysis."""

from __future__ import annotations

import os
import sqlite3
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(BASE_DIR / ".matplotlib-cache"))

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import streamlit as st


DATABASE_PATH = BASE_DIR / "cell_counts.db"
OUTPUT_DIR = BASE_DIR / "outputs"
REQUIRED_OUTPUTS = [
    "response_subject_frequencies.csv",
    "statistical_results.csv",
    "baseline_samples.csv",
    "baseline_project_counts.csv",
    "baseline_response_counts.csv",
    "baseline_sex_counts.csv",
    "final_b_cell_average.csv",
]
POPULATION_LABELS = {
    "b_cell": "B cell",
    "cd4_t_cell": "CD4 T cell",
    "cd8_t_cell": "CD8 T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}
OVERVIEW_FILTER_KEYS = [
    "overview_project",
    "overview_condition",
    "overview_treatment",
    "overview_sample_type",
    "overview_time_from_treatment_start",
]


st.set_page_config(page_title="Clinical Immune-Cell Analysis", layout="wide")


@st.cache_data
def load_overview_data(database_path: str) -> tuple[pd.DataFrame, dict[str, int]]:
    """Load dashboard metadata and per-population frequencies from SQLite."""
    with sqlite3.connect(database_path) as connection:
        query = """
            WITH totals AS (
                SELECT sample_id, SUM(cell_count) AS total_count
                FROM cell_counts
                GROUP BY sample_id
            )
            SELECT
                subjects.project_id AS project,
                subjects.condition,
                subjects.treatment,
                samples.sample_type,
                samples.time_from_treatment_start,
                samples.sample_id AS sample,
                totals.total_count,
                cell_counts.population,
                cell_counts.cell_count AS count,
                CASE WHEN totals.total_count = 0 THEN NULL
                     ELSE cell_counts.cell_count * 100.0 / totals.total_count
                END AS percentage
            FROM cell_counts
            JOIN totals ON totals.sample_id = cell_counts.sample_id
            JOIN samples ON samples.sample_id = cell_counts.sample_id
            JOIN subjects ON subjects.subject_id = samples.subject_id
            ORDER BY samples.sample_id, cell_counts.population
        """
        data = pd.read_sql_query(query, connection)
        counts = {
            table: connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            for table in ("projects", "subjects", "samples", "cell_counts")
        }
    return data, counts


@st.cache_data
def load_csv(filename: str) -> pd.DataFrame:
    """Read a generated output table and cache it for the Streamlit session."""
    return pd.read_csv(OUTPUT_DIR / filename)


def multiselect_filter(
    frame: pd.DataFrame, column: str, label: str, key: str
) -> pd.DataFrame:
    """Render an all-selected multiselect and apply it to a data frame."""
    choices = sorted(frame[column].dropna().unique().tolist())
    selected = st.multiselect(label, choices, default=choices, key=key)
    return frame[frame[column].isin(selected)]


def reset_overview_filters() -> None:
    """Return every overview filter to its default all-values selection."""
    for key in OVERVIEW_FILTER_KEYS:
        st.session_state.pop(key, None)


def render_overview() -> None:
    """Render database metrics and the filterable frequency table."""
    data, counts = load_overview_data(str(DATABASE_PATH))
    metrics = st.columns(4)
    for column, (label, table) in zip(
        metrics,
        [
            ("Projects", "projects"),
            ("Subjects", "subjects"),
            ("Samples", "samples"),
            ("Cell measurements", "cell_counts"),
        ],
    ):
        column.metric(label, f"{counts[table]:,}")

    st.subheader("Sample relative frequencies")
    st.button("Reset filters", on_click=reset_overview_filters)
    filter_columns = st.columns(5)
    filtered = data
    filter_specs = [
        ("project", "Project"),
        ("condition", "Condition"),
        ("treatment", "Treatment"),
        ("sample_type", "Sample type"),
        ("time_from_treatment_start", "Time point"),
    ]
    for container, (column, label) in zip(filter_columns, filter_specs):
        with container:
            filtered = multiselect_filter(filtered, column, label, f"overview_{column}")

    summary_column, download_column = st.columns([4, 1])
    with summary_column:
        st.caption(f"Showing {len(filtered):,} cell measurements")
    with download_column:
        st.download_button(
            "Download filtered data",
            data=filtered.to_csv(index=False).encode("utf-8"),
            file_name="filtered_relative_frequencies.csv",
            mime="text/csv",
            width="stretch",
        )

    st.dataframe(filtered, width="stretch", hide_index=True)


def render_treatment_response() -> None:
    """Render the fixed response cohort, plots, and statistical results."""
    st.info(
        "Cohort: condition = melanoma; treatment = miraclib; sample type = PBMC; "
        "response = yes or no. Each value is a subject-level mean across available time points."
    )
    values = load_csv("response_subject_frequencies.csv")
    results = load_csv("statistical_results.csv")

    response_counts = values.groupby("response")["subject"].nunique()
    left, right = st.columns(2)
    left.metric("Responder subjects", f"{response_counts.get('yes', 0):,}")
    right.metric("Non-responder subjects", f"{response_counts.get('no', 0):,}")

    population = st.selectbox(
        "Immune-cell population",
        options=list(POPULATION_LABELS),
        format_func=POPULATION_LABELS.get,
    )
    plot_data = values[values["population"] == population]
    figure, axis = plt.subplots(figsize=(7, 4.5))
    sns.boxplot(
        data=plot_data,
        x="response",
        y="mean_percentage",
        order=["no", "yes"],
        hue="response",
        palette={"no": "#7A9EAF", "yes": "#D17A57"},
        legend=False,
        ax=axis,
    )
    axis.set_title(f"{POPULATION_LABELS[population]} relative frequency")
    axis.set_xlabel("Response")
    axis.set_ylabel("Subject-level mean relative frequency (%)")
    axis.set_xticks([0, 1], labels=["Non-responder", "Responder"])
    st.pyplot(figure)
    plt.close(figure)

    st.subheader("Welch t-tests with Benjamini–Hochberg correction")
    st.dataframe(results, width="stretch", hide_index=True)
    st.download_button(
        "Download statistical results",
        data=results.to_csv(index=False).encode("utf-8"),
        file_name="statistical_results.csv",
        mime="text/csv",
    )

    st.subheader("Mean difference across populations")
    difference_data = results.copy()
    difference_data["population_label"] = difference_data["population"].map(
        POPULATION_LABELS
    )
    difference_data = difference_data.sort_values("mean_difference")
    colors = [
        "#D17A57" if value >= 0 else "#7A9EAF"
        for value in difference_data["mean_difference"]
    ]
    figure, axis = plt.subplots(figsize=(7, 3.5))
    axis.barh(
        difference_data["population_label"],
        difference_data["mean_difference"],
        color=colors,
    )
    axis.axvline(0, color="#444444", linewidth=1)
    axis.set_xlabel("Responder mean − non-responder mean (percentage points)")
    axis.set_ylabel("")
    axis.set_title("Difference in subject-level mean relative frequency")
    figure.tight_layout()
    st.pyplot(figure)
    plt.close(figure)

    with st.expander("Analysis methodology"):
        st.write(
            "Each eligible subject has repeated samples at days 0, 7, and 14. "
            "For each cell population, the dashboard uses that subject's mean relative "
            "frequency across available time points so repeated samples are not treated "
            "as independent patients."
        )
        st.write(
            "Responders and non-responders are compared with a two-sided Welch "
            "independent-samples t-test, which does not assume equal group variances. "
            "Benjamini–Hochberg correction controls the false-discovery rate across the "
            "five population tests; adjusted p-values below 0.05 are considered significant."
        )

    significant = results.loc[results["adjusted_p_value"] < 0.05, "population"].tolist()
    if significant:
        labels = ", ".join(POPULATION_LABELS.get(value, value) for value in significant)
        st.success(
            f"{labels} relative frequency is significantly associated with response in "
            "this dataset (adjusted p-value < 0.05). This association does not establish causation."
        )
    else:
        st.info("No population is significantly associated with response after FDR correction.")


def render_baseline() -> None:
    """Render the required baseline subset summaries and records."""
    baseline = load_csv("baseline_samples.csv")
    project_counts = load_csv("baseline_project_counts.csv")
    response_counts = load_csv("baseline_response_counts.csv")
    sex_counts = load_csv("baseline_sex_counts.csv")
    final_average = load_csv("final_b_cell_average.csv").loc[0, "value"]

    first, second = st.columns(2)
    first.metric("Matching baseline samples", f"{len(baseline):,}")
    second.metric("Average B-cell count", f"{final_average:,.2f}")
    st.caption(
        "The baseline cohort uses melanoma, PBMC, miraclib, and day 0. The B-cell "
        "metric separately includes responding melanoma males at day 0 across all "
        "treatments and sample types."
    )

    columns = st.columns(3)
    summaries = [
        ("Samples by project", project_counts, "project"),
        ("Subjects by response", response_counts, "response"),
        ("Subjects by sex", sex_counts, "sex"),
    ]
    for container, (title, frame, index_column) in zip(columns, summaries):
        with container:
            st.subheader(title)
            value_column = frame.columns[-1]
            st.bar_chart(frame.set_index(index_column)[value_column])
            st.dataframe(frame, width="stretch", hide_index=True)

    st.subheader("Matching baseline samples")
    st.dataframe(baseline, width="stretch", hide_index=True)
    st.download_button(
        "Download baseline records",
        data=baseline.to_csv(index=False).encode("utf-8"),
        file_name="baseline_samples.csv",
        mime="text/csv",
    )


def main() -> None:
    """Validate prerequisites and render the three dashboard sections."""
    st.title("Clinical Immune-Cell Analysis")
    st.write(
        "Explore immune-cell composition and its association with response to miraclib."
    )

    missing = []
    if not DATABASE_PATH.exists():
        missing.append(DATABASE_PATH.name)
    missing.extend(
        filename for filename in REQUIRED_OUTPUTS if not (OUTPUT_DIR / filename).exists()
    )
    if missing:
        st.error(
            "Required pipeline files are missing. Run `make pipeline` and reload this page. "
            f"Missing: {', '.join(missing)}"
        )
        st.stop()

    overview_tab, response_tab, baseline_tab = st.tabs(
        ["Data Overview", "Treatment Response", "Baseline Analysis"]
    )
    with overview_tab:
        render_overview()
    with response_tab:
        render_treatment_response()
    with baseline_tab:
        render_baseline()


if __name__ == "__main__":
    main()
