"""Run the immune-cell analyses and create all required output artifacts."""

from __future__ import annotations

import os
import sqlite3
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(BASE_DIR / ".matplotlib-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats


DATABASE_PATH = BASE_DIR / "cell_counts.db"
OUTPUT_DIR = BASE_DIR / "outputs"
POPULATION_ORDER = ["b_cell", "cd4_t_cell", "cd8_t_cell", "nk_cell", "monocyte"]
POPULATION_LABELS = {
    "b_cell": "B cell",
    "cd4_t_cell": "CD4 T cell",
    "cd8_t_cell": "CD8 T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}


def get_database_connection() -> sqlite3.Connection:
    """Open the project database with foreign-key enforcement enabled."""
    if not DATABASE_PATH.exists():
        raise FileNotFoundError("cell_counts.db is missing; run `python load_data.py` first")
    connection = sqlite3.connect(DATABASE_PATH)
    connection.execute("PRAGMA foreign_keys = ON")
    return connection


def calculate_relative_frequencies(connection: sqlite3.Connection) -> pd.DataFrame:
    """Return one count and relative-frequency row per sample population."""
    query = """
        WITH totals AS (
            SELECT sample_id, SUM(cell_count) AS total_count
            FROM cell_counts
            GROUP BY sample_id
        )
        SELECT
            counts.sample_id AS sample,
            totals.total_count,
            counts.population,
            counts.cell_count AS count,
            CASE
                WHEN totals.total_count = 0 THEN NULL
                ELSE counts.cell_count * 100.0 / totals.total_count
            END AS percentage
        FROM cell_counts AS counts
        JOIN totals ON totals.sample_id = counts.sample_id
        ORDER BY counts.sample_id, counts.population
    """
    frequencies = pd.read_sql_query(query, connection)
    expected_columns = ["sample", "total_count", "population", "count", "percentage"]
    if frequencies.columns.tolist() != expected_columns:
        raise RuntimeError("Relative-frequency output columns are incorrect")

    nonzero = frequencies[frequencies["total_count"] > 0]
    sums = nonzero.groupby("sample")["percentage"].sum()
    if not ((sums - 100).abs() < 1e-8).all():
        raise RuntimeError("Relative frequencies do not sum to 100% for every sample")
    return frequencies


def prepare_subject_level_response_data(
    connection: sqlite3.Connection, frequencies: pd.DataFrame) -> pd.DataFrame:
    """Build one response-analysis value per subject and cell population."""
    metadata_query = """
        SELECT
            samples.sample_id AS sample,
            samples.subject_id AS subject,
            subjects.response
        FROM samples
        JOIN subjects ON subjects.subject_id = samples.subject_id
        WHERE LOWER(TRIM(subjects.condition)) = 'melanoma'
          AND LOWER(TRIM(subjects.treatment)) = 'miraclib'
          AND LOWER(TRIM(samples.sample_type)) = 'pbmc'
          AND LOWER(TRIM(subjects.response)) IN ('yes', 'no')
    """
    metadata = pd.read_sql_query(metadata_query, connection)
    joined = frequencies.merge(metadata, on="sample", how="inner", validate="many_to_one")
    joined["response"] = joined["response"].str.strip().str.lower()
    # Repeated days belong to the same person, so subjects—not samples—are the units.
    subject_values = (
        joined.groupby(["subject", "response", "population"], as_index=False)[
            "percentage"
        ]
        .mean()
        .rename(columns={"percentage": "mean_percentage"})
        .sort_values(["subject", "population"])
        .reset_index(drop=True)
    )
    return subject_values


def run_response_statistics(subject_values: pd.DataFrame) -> pd.DataFrame:
    """Compare response groups with Welch tests and control the FDR at 5%."""
    rows = []
    for population in POPULATION_ORDER:
        subset = subject_values[subject_values["population"] == population]
        responders = subset.loc[subset["response"] == "yes", "mean_percentage"]
        nonresponders = subset.loc[subset["response"] == "no", "mean_percentage"]
        result = stats.ttest_ind(responders, nonresponders, equal_var=False)
        rows.append(
            {
                "population": population,
                "responder_n": len(responders),
                "non_responder_n": len(nonresponders),
                "responder_mean": responders.mean(),
                "non_responder_mean": nonresponders.mean(),
                "mean_difference": responders.mean() - nonresponders.mean(),
                "test_statistic": result.statistic,
                "p_value": result.pvalue,
            }
        )

    results = pd.DataFrame(rows)
    # Adjust the five population tests together to limit false discoveries.
    results["adjusted_p_value"] = stats.false_discovery_control(
        results["p_value"].to_numpy(), method="bh"
    )
    results["significant"] = results["adjusted_p_value"] < 0.05
    return results


def create_response_boxplots(subject_values: pd.DataFrame, output_path: Path) -> None:
    """Save a five-panel comparison of responder and non-responder values."""
    sns.set_theme(style="whitegrid", context="notebook")
    figure, axes = plt.subplots(2, 3, figsize=(13, 8), sharey=False)
    palette = {"no": "#7A9EAF", "yes": "#D17A57"}

    for axis, population in zip(axes.flat, POPULATION_ORDER):
        subset = subject_values[subject_values["population"] == population]
        sns.boxplot(
            data=subset,
            x="response",
            y="mean_percentage",
            order=["no", "yes"],
            hue="response",
            hue_order=["no", "yes"],
            palette=palette,
            legend=False,
            width=0.6,
            fliersize=2,
            ax=axis,
        )
        axis.set_title(POPULATION_LABELS[population])
        axis.set_xlabel("Response")
        axis.set_ylabel("Relative frequency (%)")
        axis.set_xticks([0, 1], labels=["Non-responder", "Responder"])
        axis.grid(axis="x", visible=False)

    axes.flat[-1].set_visible(False)
    figure.suptitle(
        "Immune-cell relative frequency by treatment response\n"
        "Subject-level means across available days 0, 7, and 14",
        fontsize=16,
        y=1.02,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(figure)


BASELINE_FROM_SQL = """
    FROM samples
    JOIN subjects ON subjects.subject_id = samples.subject_id
    WHERE LOWER(TRIM(subjects.condition)) = 'melanoma'
      AND LOWER(TRIM(samples.sample_type)) = 'pbmc'
      AND LOWER(TRIM(subjects.treatment)) = 'miraclib'
      AND samples.time_from_treatment_start = 0
"""


def query_baseline_samples(connection: sqlite3.Connection) -> pd.DataFrame:
    """Return all records in the required melanoma treatment baseline cohort."""
    query = """
        SELECT
            samples.sample_id AS sample,
            samples.subject_id AS subject,
            subjects.project_id AS project,
            subjects.condition,
            subjects.age,
            subjects.sex,
            subjects.treatment,
            subjects.response,
            samples.sample_type,
            samples.time_from_treatment_start
    """ + BASELINE_FROM_SQL + " ORDER BY samples.sample_id"
    return pd.read_sql_query(query, connection)


def query_baseline_project_counts(connection: sqlite3.Connection) -> pd.DataFrame:
    """Count baseline samples by project."""
    query = (
        "SELECT subjects.project_id AS project, COUNT(*) AS sample_count "
        + BASELINE_FROM_SQL
        + " GROUP BY subjects.project_id ORDER BY subjects.project_id"
    )
    return pd.read_sql_query(query, connection)


def query_baseline_response_counts(connection: sqlite3.Connection) -> pd.DataFrame:
    """Count distinct baseline subjects by treatment response."""
    query = (
        "SELECT LOWER(TRIM(subjects.response)) AS response, "
        "COUNT(DISTINCT subjects.subject_id) AS subject_count "
        + BASELINE_FROM_SQL
        + " GROUP BY LOWER(TRIM(subjects.response)) ORDER BY response"
    )
    return pd.read_sql_query(query, connection)


def query_baseline_sex_counts(connection: sqlite3.Connection) -> pd.DataFrame:
    """Count distinct baseline subjects by sex."""
    query = (
        "SELECT subjects.sex, COUNT(DISTINCT subjects.subject_id) AS subject_count "
        + BASELINE_FROM_SQL
        + " GROUP BY subjects.sex ORDER BY subjects.sex"
    )
    return pd.read_sql_query(query, connection)


def calculate_final_b_cell_average(
    connection: sqlite3.Connection) -> tuple[pd.DataFrame, int]:
    """Calculate the requested raw B-cell mean and qualifying subject count."""
    query = """
        SELECT
            AVG(cell_counts.cell_count) AS average_b_cell_count,
            COUNT(DISTINCT subjects.subject_id) AS qualifying_subjects
        FROM cell_counts
        JOIN samples ON samples.sample_id = cell_counts.sample_id
        JOIN subjects ON subjects.subject_id = samples.subject_id
        WHERE LOWER(TRIM(subjects.condition)) = 'melanoma'
          AND UPPER(TRIM(subjects.sex)) = 'M'
          AND LOWER(TRIM(subjects.response)) = 'yes'
          AND samples.time_from_treatment_start = 0
          AND cell_counts.population = 'b_cell'
    """
    result = pd.read_sql_query(query, connection).iloc[0]
    output = pd.DataFrame(
        {
            "metric": ["average_b_cell_count"],
            "value": [result["average_b_cell_count"]],
        }
    )
    return output, int(result["qualifying_subjects"])


def main() -> None:
    """Run all analyses and write the submission artifacts."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    with get_database_connection() as connection:
        frequencies = calculate_relative_frequencies(connection)
        subject_values = prepare_subject_level_response_data(connection, frequencies)
        statistical_results = run_response_statistics(subject_values)
        baseline_samples = query_baseline_samples(connection)
        project_counts = query_baseline_project_counts(connection)
        response_counts = query_baseline_response_counts(connection)
        sex_counts = query_baseline_sex_counts(connection)
        final_average, qualifying_subjects = calculate_final_b_cell_average(connection)

    outputs = {
        "relative_frequencies.csv": frequencies,
        "response_subject_frequencies.csv": subject_values,
        "statistical_results.csv": statistical_results,
        "baseline_samples.csv": baseline_samples,
        "baseline_project_counts.csv": project_counts,
        "baseline_response_counts.csv": response_counts,
        "baseline_sex_counts.csv": sex_counts,
        "final_b_cell_average.csv": final_average,
    }
    for filename, frame in outputs.items():
        float_format = "%.2f" if filename == "final_b_cell_average.csv" else "%.8f"
        frame.to_csv(OUTPUT_DIR / filename, index=False, float_format=float_format)

    plot_path = OUTPUT_DIR / "response_boxplots.png"
    create_response_boxplots(subject_values, plot_path)

    significant = statistical_results.loc[
        statistical_results["significant"], "population"
    ].tolist()
    eligible_subjects = subject_values["subject"].nunique()
    final_value = final_average.loc[0, "value"]
    print(f"Relative-frequency rows: {len(frequencies):,}")
    print(f"Eligible response-analysis subjects: {eligible_subjects:,}")
    print(f"Significant populations: {', '.join(significant) if significant else 'none'}")
    print(f"Baseline samples: {len(baseline_samples):,}")
    print(
        f"Average baseline B-cell count: {final_value:.2f} "
        f"({qualifying_subjects:,} qualifying subjects)"
    )
    print("Generated files:")
    for filename in [*outputs, plot_path.name]:
        print(f"  {OUTPUT_DIR / filename}")


if __name__ == "__main__":
    main()
