"""Focused integration checks for the generated database and analysis outputs."""

from __future__ import annotations

import sqlite3
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
DATABASE = ROOT / "cell_counts.db"
OUTPUTS = ROOT / "outputs"


@pytest.fixture(scope="session", autouse=True)
def completed_pipeline() -> None:
    subprocess.run([sys.executable, "load_data.py"], cwd=ROOT, check=True)
    subprocess.run([sys.executable, "analysis.py"], cwd=ROOT, check=True)


def test_database_schema_and_counts() -> None:
    with sqlite3.connect(DATABASE) as connection:
        tables = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            )
        }
        assert {"projects", "subjects", "samples", "cell_counts"} <= tables
        expected = {
            "projects": 3,
            "subjects": 3_500,
            "samples": 10_500,
            "cell_counts": 52_500,
        }
        actual = {
            table: connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            for table in expected
        }
    assert actual == expected


def test_relative_frequency_output() -> None:
    frequencies = pd.read_csv(OUTPUTS / "relative_frequencies.csv")
    assert len(frequencies) == 52_500
    assert frequencies.columns.tolist() == [
        "sample",
        "total_count",
        "population",
        "count",
        "percentage",
    ]
    representative = frequencies["sample"].drop_duplicates().iloc[[0, 100, -1]]
    sums = frequencies[frequencies["sample"].isin(representative)].groupby("sample")[
        "percentage"
    ].sum()
    assert ((sums - 100.0).abs() < 1e-5).all()


def test_baseline_subset_and_summaries() -> None:
    baseline = pd.read_csv(OUTPUTS / "baseline_samples.csv")
    projects = pd.read_csv(OUTPUTS / "baseline_project_counts.csv").set_index("project")
    responses = pd.read_csv(OUTPUTS / "baseline_response_counts.csv").set_index(
        "response"
    )
    sexes = pd.read_csv(OUTPUTS / "baseline_sex_counts.csv").set_index("sex")

    assert len(baseline) == 656
    assert projects["sample_count"].to_dict() == {"prj1": 384, "prj3": 272}
    assert responses["subject_count"].to_dict() == {"no": 325, "yes": 331}
    assert sexes["subject_count"].to_dict() == {"F": 312, "M": 344}


def test_final_b_cell_average() -> None:
    output_path = OUTPUTS / "final_b_cell_average.csv"
    result = pd.read_csv(output_path)
    assert result.loc[0, "metric"] == "average_b_cell_count"
    assert result.loc[0, "value"] == pytest.approx(10_206.15, abs=0.01)
    assert output_path.read_text().splitlines() == [
        "metric,value",
        "average_b_cell_count,10206.15",
    ]


def test_cd4_is_only_significant_population() -> None:
    results = pd.read_csv(OUTPUTS / "statistical_results.csv")
    significant = results.loc[results["significant"], "population"].tolist()
    assert significant == ["cd4_t_cell"]
